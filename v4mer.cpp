/**
 * @file v4mer.cpp
 * @brief Parallel V4 k-mer counter.
 *
 * The implementation uses the shared V4 encoding, canonicalization, hash
 * table, overflow handling, and buffered output defined in v4mer_core.hpp.
 */

#define V4MER_CORE_ONLY
#include "v4mer_core.hpp"
#undef V4MER_CORE_ONLY

#include <array>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <exception>
#include <fstream>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <sys/types.h>
#include <thread>
#include <unistd.h>

namespace parallel_v4mer {

struct InputJob {
    std::string bytes;
    size_t emit_from = 0;
    bool raw_fasta = false;
};

#ifndef V4MER_SHARD_COUNT
#define V4MER_SHARD_COUNT 64
#endif
#ifndef V4MER_UPDATE_BATCH_SIZE
#define V4MER_UPDATE_BATCH_SIZE 512
#endif
#ifndef V4MER_WINDOWS_PER_JOB
#define V4MER_WINDOWS_PER_JOB (256 * 1024)
#endif

static constexpr size_t SHARD_COUNT = V4MER_SHARD_COUNT;
static constexpr size_t MAX_SHARD_COUNT = 256;
static constexpr size_t UPDATE_BATCH_SIZE = V4MER_UPDATE_BATCH_SIZE;
static constexpr size_t MAX_PENDING_UPDATES = UPDATE_BATCH_SIZE * 4;
static constexpr size_t WINDOWS_PER_JOB = V4MER_WINDOWS_PER_JOB;
static constexpr size_t OUTPUT_BLOCK_SIZE = 1024 * 1024;
static constexpr size_t MAX_OUTPUT_FORMATTERS = 4;
static constexpr size_t OUTPUT_QUEUE_BLOCKS = 2;
static constexpr size_t DEFAULT_LOAD_PERCENT = 90;

struct MemoryPlan {
    size_t shard_count = 0;
    size_t expected_entries = 0;
    size_t capacity_total = 0;
    size_t capacity_per_shard = 0;
    size_t table_bytes = 0;
    size_t overflow_bytes = 0;
    size_t input_buffer_bytes = 0;
    size_t update_buffer_bytes = 0;
    size_t output_buffer_bytes = 0;
    size_t estimated_bytes = 0;
    size_t budget_bytes = 0;
    size_t load_percent = DEFAULT_LOAD_PERCENT;
};

static_assert(SHARD_COUNT >= 2 &&
                  (SHARD_COUNT & (SHARD_COUNT - 1)) == 0,
              "SHARD_COUNT must be a power of two");
static_assert(UPDATE_BATCH_SIZE > 0, "UPDATE_BATCH_SIZE must be positive");
static_assert(WINDOWS_PER_JOB > 0, "WINDOWS_PER_JOB must be positive");

class JobQueue {
private:
    std::deque<InputJob> jobs_;
    const size_t capacity_;
    std::mutex mutex_;
    std::condition_variable not_empty_;
    std::condition_variable not_full_;
    bool closed_ = false;
    bool cancelled_ = false;
    std::atomic<uint64_t> wait_ns_{0};

public:
    explicit JobQueue(size_t capacity) : capacity_(std::max<size_t>(capacity, 1)) {}

    bool push(InputJob job) {
        std::unique_lock<std::mutex> lock(mutex_);
        const auto wait_start = std::chrono::steady_clock::now();
        not_full_.wait(lock, [this] {
            return jobs_.size() < capacity_ || cancelled_;
        });
        wait_ns_.fetch_add(static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - wait_start).count()), std::memory_order_relaxed);
        if (cancelled_) return false;
        jobs_.push_back(std::move(job));
        not_empty_.notify_one();
        return true;
    }

    bool pop(InputJob& job) {
        std::unique_lock<std::mutex> lock(mutex_);
        const auto wait_start = std::chrono::steady_clock::now();
        not_empty_.wait(lock, [this] {
            return !jobs_.empty() || closed_ || cancelled_;
        });
        wait_ns_.fetch_add(static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - wait_start).count()), std::memory_order_relaxed);
        if (cancelled_ || jobs_.empty()) return false;
        job = std::move(jobs_.front());
        jobs_.pop_front();
        not_full_.notify_one();
        return true;
    }

    void close() {
        std::lock_guard<std::mutex> lock(mutex_);
        closed_ = true;
        not_empty_.notify_all();
    }

    void cancel() {
        std::lock_guard<std::mutex> lock(mutex_);
        cancelled_ = true;
        jobs_.clear();
        not_empty_.notify_all();
        not_full_.notify_all();
    }

    uint64_t wait_ns() const { return wait_ns_.load(std::memory_order_relaxed); }
};

class ProcessingCancelled : public std::exception {};

bool is_repetitive_job(const std::string& job) {
    static constexpr size_t SAMPLE_COUNT = 64;
    static constexpr size_t SIGNATURE_BASES = 16;
    static constexpr size_t DUPLICATE_THRESHOLD = 8;
    if (job.size() < 4096 || job.size() < SIGNATURE_BASES) return false;

    std::array<uint32_t, SAMPLE_COUNT> signatures{};
    size_t duplicates = 0;
    const size_t last_start = job.size() - SIGNATURE_BASES;
    for (size_t sample = 0; sample < SAMPLE_COUNT; ++sample) {
        const size_t start = (sample * last_start) / (SAMPLE_COUNT - 1);
        uint32_t signature = 0;
        for (size_t base = 0; base < SIGNATURE_BASES; ++base) {
            signature = (signature << 2) |
                        static_cast<unsigned char>(job[start + base]);
        }
        for (size_t previous = 0; previous < sample; ++previous) {
            if (signatures[previous] == signature) {
                ++duplicates;
                break;
            }
        }
        if (duplicates >= DUPLICATE_THRESHOLD) return true;
        signatures[sample] = signature;
    }
    return false;
}

class RawOutputFile {
private:
    FILE* file_ = nullptr;

public:
    explicit RawOutputFile(const std::string& filename) {
        file_ = fopen(filename.c_str(), "wb");
        if (!file_) {
            throw std::runtime_error("Cannot open output file: " + filename +
                                     ": " + std::strerror(errno));
        }
    }

    ~RawOutputFile() {
        if (file_) fclose(file_);
    }

    RawOutputFile(const RawOutputFile&) = delete;
    RawOutputFile& operator=(const RawOutputFile&) = delete;

    void write(const std::string& block) {
        if (fwrite(block.data(), 1, block.size(), file_) != block.size()) {
            throw std::runtime_error("Error writing output file");
        }
    }

    void close() {
        FILE* file = file_;
        file_ = nullptr;
        if (fclose(file) != 0) {
            throw std::runtime_error("Error closing output file");
        }
    }
};

class SequenceChunker {
private:
    JobQueue& queue_;
    const size_t k_;
    std::string run_;

    void enqueue(std::string job) {
        if (!queue_.push(InputJob{std::move(job), 0, false})) throw ProcessingCancelled();
    }

public:
    SequenceChunker(JobQueue& queue, size_t k) : queue_(queue), k_(k) {
        run_.reserve(WINDOWS_PER_JOB + k - 1);
    }

    void add_base(uint8_t encoded) {
        run_.push_back(static_cast<char>(encoded));
        const size_t job_size = WINDOWS_PER_JOB + k_ - 1;
        if (run_.size() == job_size) {
            std::string overlap;
            if (k_ > 1) {
                overlap.assign(run_.end() - static_cast<std::ptrdiff_t>(k_ - 1),
                               run_.end());
            }
            enqueue(std::move(run_));
            run_ = std::move(overlap);
            run_.reserve(job_size);
        }
    }

    void break_run() {
        if (run_.size() >= k_) enqueue(std::move(run_));
        run_.clear();
        run_.reserve(WINDOWS_PER_JOB + k_ - 1);
    }

    void finish() { break_run(); }
};

class ParallelSequenceReader {
private:
    std::string filename_;
    CompressionType compression_;
    FileFormat format_;

    void produce_fasta(BufferedReader& reader, SequenceChunker& chunker) {
        bool in_header = false;
        bool at_line_start = true;
        int c;
        while ((c = reader.get()) != EOF) {
            if (in_header) {
                if (c == '\n') {
                    in_header = false;
                    at_line_start = true;
                }
                continue;
            }
            if (c == '\n') {
                at_line_start = true;
                continue;
            }
            if (c == '\r' || c == ' ' || c == '\t') continue;
            if (at_line_start && (c == '>' || c == ';')) {
                in_header = true;
                if (c == '>') chunker.break_run();
                continue;
            }

            at_line_start = false;
            const int8_t code = BASE_ENCODING[static_cast<unsigned char>(c)];
            if (code >= 0) {
                chunker.add_base(static_cast<uint8_t>(code));
            } else {
                chunker.break_run();
            }
        }
        chunker.finish();
    }

    void produce_plain_fasta(size_t k, JobQueue& queue) {
        std::ifstream input(filename_);
        if (!input) throw std::runtime_error("Cannot open file: " + filename_);
        const size_t overlap_limit = k - 1 + 1024;
        std::string pending;
        std::string overlap;
        pending.reserve(WINDOWS_PER_JOB + overlap_limit);
        std::string line;
        while (std::getline(input, line)) {
            pending.append(line);
            pending.push_back(10);
            if (pending.size() < WINDOWS_PER_JOB) continue;
            InputJob job;
            job.bytes = overlap + pending;
            job.emit_from = overlap.size();
            job.raw_fasta = true;
            if (!queue.push(std::move(job))) throw ProcessingCancelled();
            size_t keep_start = pending.size() > overlap_limit
                                    ? pending.size() - overlap_limit
                                    : 0;
            while (keep_start > 0 && pending[keep_start - 1] != 10) --keep_start;
            overlap.assign(pending, keep_start, std::string::npos);
            pending.clear();
        }
        if (input.bad()) throw std::runtime_error("Error reading file: " + filename_);
        if (!pending.empty()) {
            InputJob job;
            job.bytes = overlap + pending;
            job.emit_from = overlap.size();
            job.raw_fasta = true;
            if (!queue.push(std::move(job))) throw ProcessingCancelled();
        }
    }

    void produce_fastq(BufferedReader& reader, SequenceChunker& chunker) {
        while (!reader.eof()) {
            int c = reader.peek();
            while (c == '\n' || c == '\r') {
                reader.get();
                c = reader.peek();
            }
            if (c == EOF) break;
            if (c != '@') {
                throw std::runtime_error("Malformed FASTQ: expected '@' header");
            }

            reader.skip_line();
            chunker.break_run();
            size_t sequence_length = 0;

            while (true) {
                c = reader.peek();
                if (c == EOF) {
                    throw std::runtime_error("Malformed FASTQ: missing '+' line");
                }
                if (c == '+') {
                    reader.skip_line();
                    chunker.break_run();
                    break;
                }

                while ((c = reader.get()) != '\n' && c != EOF) {
                    if (c == '\r') continue;
                    if (sequence_length == std::numeric_limits<size_t>::max()) {
                        throw std::length_error("FASTQ sequence is too long");
                    }
                    ++sequence_length;
                    const int8_t code =
                        BASE_ENCODING[static_cast<unsigned char>(c)];
                    if (code >= 0) {
                        chunker.add_base(static_cast<uint8_t>(code));
                    } else {
                        chunker.break_run();
                    }
                }
                if (c == EOF) {
                    throw std::runtime_error("Malformed FASTQ: missing '+' line");
                }
            }

            size_t quality_length = 0;
            while (quality_length < sequence_length) {
                c = reader.get();
                if (c == EOF) {
                    throw std::runtime_error(
                        "Malformed FASTQ: truncated quality data");
                }
                if (c != '\n' && c != '\r') ++quality_length;
            }

            c = reader.peek();
            while (c == '\r') {
                reader.get();
                c = reader.peek();
            }
            if (c == '\n') {
                reader.get();
            } else if (c != EOF) {
                throw std::runtime_error(
                    "Malformed FASTQ: quality length exceeds sequence length");
            }
        }
        chunker.finish();
    }

public:
    explicit ParallelSequenceReader(const std::string& filename)
        : filename_(filename), compression_(detect_compression(filename)) {
        BufferedReader reader(filename_, compression_);
        format_ = detect_format(reader);
    }

    FileFormat format() const { return format_; }
    CompressionType compression() const { return compression_; }

    void produce(size_t k, JobQueue& queue) {
        if (compression_ == CompressionType::NONE && format_ == FileFormat::FASTA) {
            produce_plain_fasta(k, queue);
            return;
        }
        BufferedReader reader(filename_, compression_);
        SequenceChunker chunker(queue, k);
        if (format_ == FileFormat::FASTA) {
            produce_fasta(reader, chunker);
        } else {
            produce_fastq(reader, chunker);
        }
    }
};

template<size_t WORDS>
class ShardedV4Table {
public:
    struct PendingUpdate {
        std::array<uint64_t, WORDS> kmer{};
        Transform transform = TRANSFORM_I;
        uint32_t count = 1;
        size_t hash = 0;
    };

    struct __attribute__((packed)) OutputRecord {
        // kmer is already in the rendered lexicographic orientation.
        std::array<uint64_t, WORDS> kmer{};
        uint64_t count = 0;
    };

    struct WorkerBuffers {
        std::vector<std::vector<PendingUpdate>> updates;
        explicit WorkerBuffers(size_t shard_count) : updates(shard_count) {}
        size_t pending_count = 0;
    };

private:
    using Table = CompactHashTable<WORDS>;

    struct alignas(64) Shard {
        std::mutex mutex;
        std::unique_ptr<Table> table;
    };

    std::array<Shard, MAX_SHARD_COUNT> shards_;
    size_t shard_count_;
    unsigned shard_bits_;
    std::atomic<uint64_t>* mutex_wait_ns_ = nullptr;

    static constexpr unsigned compute_shard_bits(size_t count) {
        unsigned bits = 0;
        while (count > 1) {
            ++bits;
            count >>= 1;
        }
        return bits;
    }

    size_t shard_index(size_t hash) const {
        return static_cast<size_t>(hash >> (64 - shard_bits_));
    }

    static bool keys_equal(const uint64_t* left, const uint64_t* right) {
        for (size_t word = 0; word < WORDS; ++word) {
            if (left[word] != right[word]) return false;
        }
        return true;
    }

    static std::array<uint64_t, WORDS> make_output_kmer(
        const uint64_t* canonical, size_t k, OutputTransform transform) {
        std::array<uint64_t, WORDS> rendered{};
        for (size_t position = 0; position < k; ++position) {
            const size_t source =
                transform == OutputTransform::REVERSE ? k - 1 - position : position;
            uint8_t base = encoded_base(canonical, source);
            if (transform == OutputTransform::COMPLEMENT) base ^= 3;
            const size_t word = position >> 5;
            const size_t offset = (position & 31) << 1;
            rendered[word] |= static_cast<uint64_t>(base) << offset;
        }
        return rendered;
    }

    template<bool BULK>
    void flush(size_t index, WorkerBuffers& buffers) {
        auto& updates = buffers.updates[index];
        if (updates.empty()) return;
        Shard& shard = shards_[index];
        const auto wait_start = std::chrono::steady_clock::now();
        std::unique_lock<std::mutex> lock(shard.mutex);
        if (mutex_wait_ns_ != nullptr) {
            mutex_wait_ns_->fetch_add(static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - wait_start).count()), std::memory_order_relaxed);
        }
        for (const PendingUpdate& update : updates) {
            if constexpr (!BULK) {
                shard.table->insert_or_add_hashed(update.kmer.data(),
                                                  update.transform, 1,
                                                  update.hash);
            } else {
                shard.table->insert_or_add_hashed(update.kmer.data(),
                                                  update.transform,
                                                  update.count, update.hash);
            }
        }
        buffers.pending_count -= updates.size();
        updates.clear();
    }

    static void append_output_line(std::string& block,
                                   const uint64_t* kmer_data, size_t k,
                                   OutputTransform transform, uint64_t count) {
        for (size_t i = 0; i < k; ++i) {
            const size_t source = transform == OutputTransform::REVERSE ? k - 1 - i : i;
            uint8_t base = encoded_base(kmer_data, source);
            if (transform == OutputTransform::COMPLEMENT) base ^= 3;
            block.push_back(BASE_DECODING[base]);
        }
        block.push_back(9);
        char count_buffer[32];
        const auto conversion = std::to_chars(count_buffer, count_buffer + sizeof(count_buffer), count);
        if (conversion.ec != std::errc()) throw std::runtime_error("Cannot format k-mer count");
        block.append(count_buffer, conversion.ptr);
        block.push_back(10);
    }

    void format_shard(size_t index, size_t k,
                      std::vector<std::vector<std::string>>& output_blocks,
                      std::mutex& output_mutex,
                      std::atomic<size_t>& output_lines) const {
        const Table& table = *shards_[index].table;
        std::vector<std::string> lines;
        lines.reserve(table.size() * 2);
        size_t local_lines = 0;
        for (const auto& entry : table) {
            CompactKmer<WORDS> canonical;
            for (size_t word = 0; word < WORDS; ++word) canonical.data()[word] = entry.kmer[word];
            const uint64_t count_I = table.get_count(entry, TRANSFORM_I);
            const uint64_t count_R = table.get_count(entry, TRANSFORM_R);
            const uint64_t count_C = table.get_count(entry, TRANSFORM_C);
            const uint64_t count_RC = table.get_count(entry, TRANSFORM_RC);
            const uint64_t count_pair1 = count_I + count_RC;
            const uint64_t count_pair2 = count_R + count_C;
            if (count_pair1 > 0) {
                std::string line;
                append_output_line(line, canonical.data(), k, OutputTransform::IDENTITY, count_pair1);
                lines.push_back(std::move(line));
                ++local_lines;
            }
            if (count_pair2 > 0) {
                const OutputTransform pair_transform = reverse_is_canonical_pair_member(canonical.data(), k) ? OutputTransform::REVERSE : OutputTransform::COMPLEMENT;
                std::string line;
                append_output_line(line, canonical.data(), k, pair_transform, count_pair2);
                lines.push_back(std::move(line));
                ++local_lines;
            }
        }
        std::sort(lines.begin(), lines.end());
        std::string block;
        block.reserve(OUTPUT_BLOCK_SIZE + 256);
        for (std::string& line : lines) {
            if (block.size() + line.size() > OUTPUT_BLOCK_SIZE && !block.empty()) {
                std::lock_guard<std::mutex> lock(output_mutex);
                output_blocks[index].push_back(std::move(block));
                block.clear();
                block.reserve(OUTPUT_BLOCK_SIZE + 256);
            }
            block.append(std::move(line));
        }
        if (!block.empty()) {
            std::lock_guard<std::mutex> lock(output_mutex);
            output_blocks[index].push_back(std::move(block));
        }
        output_lines.fetch_add(local_lines, std::memory_order_relaxed);
    }

    template<typename EmitBlock>
    void format_shard_compact(size_t index, size_t k,
                              size_t& output_lines,
                              double& formatting_seconds,
                              EmitBlock&& emit_block) const {
        const auto format_start = std::chrono::steady_clock::now();
        const Table& table = *shards_[index].table;
        std::vector<OutputRecord> records;
        records.reserve(table.size() * 2);
        for (const auto& entry : table) {
            OutputRecord first;
            for (size_t word = 0; word < WORDS; ++word) {
                first.kmer[word] = entry.kmer[word];
            }
            first.count = table.get_count(entry, TRANSFORM_I) +
                          table.get_count(entry, TRANSFORM_RC);
            if (first.count > 0) records.push_back(first);

            const OutputTransform pair_transform =
                reverse_is_canonical_pair_member(first.kmer.data(), k)
                    ? OutputTransform::REVERSE
                    : OutputTransform::COMPLEMENT;
            OutputRecord second;
            second.kmer = make_output_kmer(first.kmer.data(), k, pair_transform);
            second.count = table.get_count(entry, TRANSFORM_R) +
                           table.get_count(entry, TRANSFORM_C);
            if (second.count > 0) records.push_back(second);
        }

        std::sort(records.begin(), records.end(),
                  [](const OutputRecord& left, const OutputRecord& right) {
            for (size_t word = 0; word < WORDS; ++word) {
                const uint64_t differing_bits =
                    left.kmer[word] ^ right.kmer[word];
                if (differing_bits != 0) {
                    const unsigned bit =
                        static_cast<unsigned>(__builtin_ctzll(differing_bits)) & ~1U;
                    const uint8_t lhs =
                        static_cast<uint8_t>((left.kmer[word] >> bit) & 3);
                    const uint8_t rhs =
                        static_cast<uint8_t>((right.kmer[word] >> bit) & 3);
                    return lhs < rhs;
                }
            }
            return left.count < right.count;
        });

        std::string block;
        block.reserve(OUTPUT_BLOCK_SIZE + 256);
        size_t block_lines = 0;
        for (const OutputRecord& record : records) {
            append_output_line(block, record.kmer.data(), k,
                               OutputTransform::IDENTITY, record.count);
            ++block_lines;
            ++output_lines;
            if (block.size() >= OUTPUT_BLOCK_SIZE) {
                emit_block(std::move(block), block_lines, false);
                block = std::string();
                block.reserve(OUTPUT_BLOCK_SIZE + 256);
                block_lines = 0;
            }
        }
        emit_block(std::move(block), block_lines, true);
        formatting_seconds += std::chrono::duration<double>(
            std::chrono::steady_clock::now() - format_start).count();
    }
public:
    explicit ShardedV4Table(size_t expected_entries, const MemoryPlan& plan,
                            std::atomic<uint64_t>* mutex_wait_ns)
        : shard_count_(plan.shard_count),
          shard_bits_(compute_shard_bits(plan.shard_count)),
          mutex_wait_ns_(mutex_wait_ns) {
        if (shard_count_ < 2 || shard_count_ > MAX_SHARD_COUNT ||
            (shard_count_ & (shard_count_ - 1)) != 0) {
            throw std::invalid_argument("Shard count must be a power of two");
        }
        const size_t per_shard = expected_entries / shard_count_ +
                                 (expected_entries % shard_count_ != 0 ? 1 : 0);
        for (size_t index = 0; index < shard_count_; ++index) {
            shards_[index].table =
                std::make_unique<Table>(per_shard, plan.capacity_per_shard,
                                         plan.load_percent);
        }
    }

    template<bool COALESCE_REPEATS>
    void add(const uint64_t* kmer_data, Transform transform,
             WorkerBuffers& buffers) {
        const size_t hash = Table::compute_hash(kmer_data);
        const size_t index = shard_index(hash);
        std::vector<PendingUpdate>& updates = buffers.updates[index];
        if (updates.capacity() == 0) {
            updates.reserve(std::min<size_t>(UPDATE_BATCH_SIZE, 128));
        }

        // Repetitive sequence usually produces a contiguous run of the same
        // canonical key in a shard. Fold its transform channels in the worker
        // buffer so millions of observations become a handful of locked probes.
        if constexpr (COALESCE_REPEATS) {
            if (!updates.empty() &&
                keys_equal(updates.back().kmer.data(), kmer_data)) {
                bool saturated = false;
                for (auto update = updates.rbegin(); update != updates.rend();
                     ++update) {
                    if (!keys_equal(update->kmer.data(), kmer_data)) break;
                    if (update->transform == transform) {
                        if (update->count ==
                            std::numeric_limits<uint32_t>::max()) {
                            // Preserve an unbounded final count without
                            // enlarging every pending update.
                            saturated = true;
                            break;
                        }
                        ++update->count;
                        return;
                    }
                }
                if (saturated) flush<true>(index, buffers);
            }
        }

        PendingUpdate update;
        for (size_t word = 0; word < WORDS; ++word) {
            update.kmer[word] = kmer_data[word];
        }
        update.transform = transform;
        update.hash = hash;
        if (buffers.pending_count >= MAX_PENDING_UPDATES) {
            size_t largest = 0;
            for (size_t candidate = 1; candidate < buffers.updates.size(); ++candidate) {
                if (buffers.updates[candidate].size() > buffers.updates[largest].size()) {
                    largest = candidate;
                }
            }
            flush<COALESCE_REPEATS>(largest, buffers);
        }
        updates.push_back(update);
        ++buffers.pending_count;
        if (updates.size() >= UPDATE_BATCH_SIZE) {
            flush<COALESCE_REPEATS>(index, buffers);
        }
    }

    template<bool BULK>
    void flush_all(WorkerBuffers& buffers) {
        for (size_t index = 0; index < shard_count_; ++index) {
            flush<BULK>(index, buffers);
        }
    }

    size_t size() const {
        size_t result = 0;
        for (size_t index = 0; index < shard_count_; ++index) {
            result += shards_[index].table->size();
        }
        return result;
    }

    size_t capacity() const {
        size_t result = 0;
        for (size_t index = 0; index < shard_count_; ++index) {
            result += shards_[index].table->capacity();
        }
        return result;
    }

    size_t overflow_size() const {
        size_t result = 0;
        for (size_t index = 0; index < shard_count_; ++index) {
            result += shards_[index].table->overflow_size();
        }
        return result;
    }

    size_t write_output(const std::string& output_file, size_t k,
                        size_t thread_count, double& formatting_seconds,
                        double& writing_seconds) const {
        const size_t formatter_count =
            std::min({thread_count, shard_count_, MAX_OUTPUT_FORMATTERS});

        struct OutputBlock {
            std::string data;
            size_t lines = 0;
            bool last = false;
        };
        using BlockKey = std::pair<size_t, size_t>;
        std::map<BlockKey, OutputBlock> pending_blocks;
        const size_t queue_limit = std::max<size_t>(OUTPUT_QUEUE_BLOCKS, 2);
        std::vector<size_t> shard_lines(shard_count_, 0);
        std::vector<double> shard_formatting(shard_count_, 0.0);
        std::mutex output_mutex;
        std::condition_variable output_cv;
        size_t next_claim = 0;
        size_t next_output_shard = 0;
        size_t queued_blocks = 0;
        std::exception_ptr formatting_error;
        std::vector<std::thread> formatters;

        auto formatter = [&] {
            try {
                while (true) {
                    size_t index = 0;
                    {
                        std::lock_guard<std::mutex> lock(output_mutex);
                        if (formatting_error || next_claim >= shard_count_) return;
                        index = next_claim++;
                    }

                    size_t next_block = 0;
                    auto emit_block = [&](std::string block, size_t lines,
                                          bool last) {
                        std::unique_lock<std::mutex> lock(output_mutex);
                        output_cv.wait(lock, [&] {
                            return formatting_error || queued_blocks < queue_limit ||
                                   index == next_output_shard;
                        });
                        if (formatting_error) {
                            throw std::runtime_error("Output formatting cancelled");
                        }
                        pending_blocks.emplace(
                            BlockKey{index, next_block++},
                            OutputBlock{std::move(block), lines, last});
                        ++queued_blocks;
                        lock.unlock();
                        output_cv.notify_all();
                    };
                    format_shard_compact(index, k, shard_lines[index],
                                         shard_formatting[index], emit_block);
                    output_cv.notify_all();
                }
            } catch (...) {
                {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    if (!formatting_error) formatting_error = std::current_exception();
                }
                output_cv.notify_all();
            }
        };

        formatters.reserve(formatter_count);
        for (size_t i = 0; i < formatter_count; ++i) {
            formatters.emplace_back(formatter);
        }

        RawOutputFile output(output_file);
        size_t total_lines = 0;
        for (size_t index = 0; index < shard_count_; ++index) {
            size_t block_index = 0;
            bool last = false;
            while (!last) {
                OutputBlock block;
                {
                    std::unique_lock<std::mutex> lock(output_mutex);
                    output_cv.wait(lock, [&] {
                        return formatting_error ||
                               pending_blocks.find(BlockKey{index, block_index}) !=
                                   pending_blocks.end();
                    });
                    if (formatting_error) break;
                    auto found = pending_blocks.find(BlockKey{index, block_index});
                    block = std::move(found->second);
                    pending_blocks.erase(found);
                    --queued_blocks;
                }
                output_cv.notify_all();

                total_lines += block.lines;
                const auto write_start = std::chrono::steady_clock::now();
                output.write(block.data);
                writing_seconds += std::chrono::duration<double>(
                    std::chrono::steady_clock::now() - write_start).count();
                last = block.last;
                ++block_index;
                if (last) {
                    std::lock_guard<std::mutex> lock(output_mutex);
                    next_output_shard = index + 1;
                    output_cv.notify_all();
                }
            }
            if (formatting_error) break;
        }

        for (std::thread& thread : formatters) thread.join();
        output.close();
        if (formatting_error) std::rethrow_exception(formatting_error);
        for (double seconds : shard_formatting) formatting_seconds += seconds;
        return total_lines;
    }
#if 0
            const auto write_start = std::chrono::steady_clock::now();
        const size_t formatter_count =
            std::min({thread_count, shard_count_, MAX_OUTPUT_FORMATTERS});
        std::vector<std::vector<std::string>> output_blocks(shard_count_);
        std::mutex output_mutex;
        std::atomic<size_t> next_shard{0};
        std::atomic<size_t> output_lines{0};
        std::mutex error_mutex;
        std::exception_ptr formatter_error;

        auto formatter = [&] {
            try {
                while (true) {
                    const size_t index =
                        next_shard.fetch_add(1, std::memory_order_relaxed);
                    if (index >= shard_count_) break;
                    format_shard(index, k, output_blocks, output_mutex,
                                 output_lines);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!formatter_error) formatter_error = std::current_exception();
            }
        };

        std::vector<std::thread> formatters;
        formatters.reserve(formatter_count);
        for (size_t i = 0; i < formatter_count; ++i) {
            formatters.emplace_back(formatter);
        }
        for (std::thread& thread : formatters) thread.join();
        if (formatter_error) std::rethrow_exception(formatter_error);

        RawOutputFile output(output_file);
        for (size_t index = 0; index < shard_count_; ++index) {
            for (const std::string& block : output_blocks[index]) {
                output.write(block);
            }
        }
        output.close();
        return output_lines.load(std::memory_order_relaxed);
#endif
};

template<size_t WORDS, bool COALESCE_REPEATS>
void process_job(const std::string& job, CompactRollingWindow<WORDS>& window,
                 ShardedV4Table<WORDS>& table,
                 typename ShardedV4Table<WORDS>::WorkerBuffers& buffers) {
    window.reset();
    for (unsigned char encoded : job) {
        if (window.add_base(encoded)) {
            auto canonical = window.get_canonical();
            table.template add<COALESCE_REPEATS>(
                canonical.kmer.data(), canonical.transform, buffers);
        }
    }
}

template<size_t WORDS, bool COALESCE_REPEATS>
void process_raw_job(const InputJob& job,
                     CompactRollingWindow<WORDS>& window,
                     ShardedV4Table<WORDS>& table,
                     typename ShardedV4Table<WORDS>::WorkerBuffers& buffers) {
    window.reset();
    bool in_header = false;
    bool at_line_start = true;
    for (size_t index = 0; index < job.bytes.size(); ++index) {
        const unsigned char c = static_cast<unsigned char>(job.bytes[index]);
        if (in_header) {
            if (c == 10) { in_header = false; at_line_start = true; }
            continue;
        }
        if (c == 10) { at_line_start = true; continue; }
        if (c == 13 || c == 32 || c == 9) continue;
        if (at_line_start && (c == 62 || c == 59)) {
            in_header = true;
            if (c == 62) window.reset();
            continue;
        }
        at_line_start = false;
        const int8_t code = BASE_ENCODING[c];
        if (code < 0) { window.reset(); continue; }
        const bool ready = window.add_base(static_cast<uint8_t>(code));
        if (index >= job.emit_from && ready) {
            auto canonical = window.get_canonical();
            table.template add<COALESCE_REPEATS>(
                canonical.kmer.data(), canonical.transform, buffers);
        }
    }
}

size_t choose_shard_count_by_estimate(size_t estimated_entries) {
    size_t target = 8;
    while (target < MAX_SHARD_COUNT &&
           target <= std::numeric_limits<size_t>::max() / 4096 &&
           target * 4096 < estimated_entries / 2) {
        target <<= 1;
    }
    return target;
}

size_t choose_shard_count(size_t thread_count) {
    (void)thread_count;
    // A stable shard topology is part of the byte-for-byte output contract.
    // It also keeps hash routing comparable across thread counts.
    return 64;
}

size_t automatic_memory_budget() {
    const long pages = sysconf(_SC_PHYS_PAGES);
    const long page_size = sysconf(_SC_PAGE_SIZE);
    if (pages > 0 && page_size > 0) {
        const size_t physical = static_cast<size_t>(pages) *
                                static_cast<size_t>(page_size);
        return physical / 2;
    }
    return static_cast<size_t>(4ULL * 1024ULL * 1024ULL * 1024ULL);
}

MemoryPlan make_memory_plan(size_t expected_entries, size_t thread_count,
                            size_t entry_size, size_t pending_size,
                            size_t output_record_size,
                            size_t memory_budget_mb) {
    MemoryPlan plan;
    plan.shard_count = choose_shard_count(thread_count);
    plan.expected_entries = expected_entries;
    const size_t required = expected_entries == 0
                                ? plan.shard_count * 64
                                : (expected_entries * 100 + plan.load_percent - 1) / plan.load_percent;
    plan.capacity_total = std::max<size_t>(required, plan.shard_count * 64);
    plan.capacity_per_shard = (plan.capacity_total + plan.shard_count - 1) /
                              plan.shard_count;
    plan.capacity_total = plan.capacity_per_shard * plan.shard_count;
    plan.table_bytes = plan.capacity_total * entry_size;
    plan.overflow_bytes = std::max<size_t>(sizeof(uint64_t) * 8,
                                           expected_entries / 20);
    const size_t queue_capacity = std::max<size_t>(thread_count * 4, 8);
    plan.input_buffer_bytes = queue_capacity * WINDOWS_PER_JOB;
    plan.update_buffer_bytes = thread_count * plan.shard_count *
                               MAX_PENDING_UPDATES * pending_size;
    const size_t formatter_count =
        std::min({thread_count, plan.shard_count, MAX_OUTPUT_FORMATTERS});
    const size_t records_per_shard =
        ((expected_entries + plan.shard_count - 1) / plan.shard_count) * 2;
    plan.output_buffer_bytes =
        formatter_count * records_per_shard * output_record_size +
        OUTPUT_QUEUE_BLOCKS * OUTPUT_BLOCK_SIZE;
    plan.estimated_bytes = plan.table_bytes + plan.overflow_bytes +
                           plan.input_buffer_bytes + plan.update_buffer_bytes +
                           plan.output_buffer_bytes;
    plan.budget_bytes = memory_budget_mb == 0
                            ? automatic_memory_budget()
                            : memory_budget_mb * 1024ULL * 1024ULL;
    return plan;
}

template<size_t WORDS>
int run_parallel(const std::string& input_file, size_t k,
                 const std::string& output_file, size_t file_size,
                 size_t thread_count, size_t memory_budget_mb) {
    ParallelSequenceReader reader(input_file);
    const char* format = reader.format() == FileFormat::FASTA ? "FASTA" : "FASTQ";
    const char* compression =
        reader.compression() == CompressionType::GZIP ? " (gzip)" : "";
    const size_t estimated = estimate_unique_kmers(
        file_size, k, reader.format(), reader.compression());

    std::cerr << "Format: " << format << compression << '\n';
    std::cerr << "Workers: " << thread_count << '\n';
    const MemoryPlan plan = make_memory_plan(
        estimated, thread_count, sizeof(typename CompactHashTable<WORDS>::Entry),
        sizeof(typename ShardedV4Table<WORDS>::PendingUpdate),
        sizeof(typename ShardedV4Table<WORDS>::OutputRecord),
        memory_budget_mb);
    if (plan.estimated_bytes > plan.budget_bytes) {
        throw std::runtime_error(
            "Memory budget exceeded: estimated " +
            std::to_string(plan.estimated_bytes / (1024 * 1024)) +
            " MiB, limit " + std::to_string(plan.budget_bytes / (1024 * 1024)) +
            " MiB; shards=" + std::to_string(plan.shard_count) +
            ", capacity=" + std::to_string(plan.capacity_total) +
            ". Reduce threads or use a smaller input.");
    }
    const size_t shard_count = plan.shard_count;
    std::cerr << "Memory plan: estimated " << plan.estimated_bytes / (1024 * 1024)
              << " MiB, budget " << plan.budget_bytes / (1024 * 1024)
              << " MiB, capacity " << plan.capacity_total << '\n';
    std::cerr << "Estimated table memory: " << plan.table_bytes / (1024 * 1024)
              << " MiB\n";
    std::cerr << "Output buffer memory: " << plan.output_buffer_bytes / (1024 * 1024)
              << " MiB\n";
    std::cerr << "Hash shards: " << shard_count << '\n';
    std::cerr << "Estimated unique k-mers: " << estimated << '\n';

    std::atomic<uint64_t> shard_mutex_wait_ns{0};
    ShardedV4Table<WORDS> table(estimated, plan, &shard_mutex_wait_ns);
    const size_t queue_capacity = std::max<size_t>(thread_count * 4, 8);
    JobQueue queue(queue_capacity);
    std::mutex error_mutex;
    std::exception_ptr worker_error;
    std::atomic<size_t> coalesced_jobs{0};

    auto worker = [&] {
        try {
            CompactRollingWindow<WORDS> window(k);
            typename ShardedV4Table<WORDS>::WorkerBuffers buffers(shard_count);
            bool have_buffer_mode = false;
            bool bulk_buffer_mode = false;
            InputJob job;
            while (queue.pop(job)) {
                const bool repetitive = is_repetitive_job(job.bytes);
                if (have_buffer_mode && repetitive != bulk_buffer_mode) {
                    if (bulk_buffer_mode) {
                        table.template flush_all<true>(buffers);
                    } else {
                        table.template flush_all<false>(buffers);
                    }
                }
                have_buffer_mode = true;
                bulk_buffer_mode = repetitive;

                if (repetitive) {
                    coalesced_jobs.fetch_add(1, std::memory_order_relaxed);
                    if (job.raw_fasta) {
                        process_raw_job<WORDS, true>(job, window, table, buffers);
                    } else {
                        process_job<WORDS, true>(job.bytes, window, table, buffers);
                    }
                } else if (job.raw_fasta) {
                    process_raw_job<WORDS, false>(job, window, table, buffers);
                } else {
                    process_job<WORDS, false>(job.bytes, window, table, buffers);
                }
            }
            if (have_buffer_mode) {
                if (bulk_buffer_mode) {
                    table.template flush_all<true>(buffers);
                } else {
                    table.template flush_all<false>(buffers);
                }
            }
        } catch (...) {
            {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!worker_error) worker_error = std::current_exception();
            }
            queue.cancel();
        }
    };

    const auto counting_start = std::chrono::steady_clock::now();
    std::vector<std::thread> workers;
    workers.reserve(thread_count);
    for (size_t i = 0; i < thread_count; ++i) workers.emplace_back(worker);

    std::exception_ptr producer_error;
    try {
        reader.produce(k, queue);
        queue.close();
    } catch (const ProcessingCancelled&) {
        queue.cancel();
    } catch (...) {
        producer_error = std::current_exception();
        queue.cancel();
    }

    for (std::thread& thread : workers) thread.join();
    if (worker_error) std::rethrow_exception(worker_error);
    if (producer_error) std::rethrow_exception(producer_error);

    const auto counting_end = std::chrono::steady_clock::now();
    const double counting_seconds =
        std::chrono::duration<double>(counting_end - counting_start).count();
    std::cerr << "Queue wait time: " << queue.wait_ns() / 1e9 << " s\n";
    std::cerr << "Shard mutex wait time: " << shard_mutex_wait_ns.load() / 1e9 << " s\n";

    std::cerr << "Found " << table.size() << " distinct V4 classes\n";
    std::cerr << "Hash table capacity: " << table.capacity() << " entries\n";
    std::cerr << "Hash table load factor: "
              << (100.0 * static_cast<double>(table.size()) /
                  static_cast<double>(table.capacity()))
              << "%\n";
    if (table.overflow_size() > 0) {
        std::cerr << "Overflow entries: " << table.overflow_size() << '\n';
    }
    std::cerr << "Coalesced low-complexity jobs: "
              << coalesced_jobs.load(std::memory_order_relaxed) << '\n';
    std::cerr << "Parallel counting time: " << counting_seconds << " s\n";
    std::cerr << "Writing output...\n";
    std::cerr << "Output sort: packed-2bit key\n";
    std::cerr << "Output formatting: direct block append\n";

    const auto output_start = std::chrono::steady_clock::now();
    double formatting_seconds = 0.0;
    double writing_seconds = 0.0;
    const size_t output_lines = table.write_output(
        output_file, k, thread_count, formatting_seconds, writing_seconds);
    const double output_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - output_start).count();
    std::cerr << "Wrote " << output_lines << " output lines\n";
    std::cerr << "Output time: " << output_seconds << " s\n";
    std::cerr << "Formatting time: " << formatting_seconds << " s\n";
    std::cerr << "Writing time: " << writing_seconds << " s\n";
    return 0;
}

bool parse_size(const char* text, size_t& value) {
    const char* end = text + std::strlen(text);
    const auto parsed = std::from_chars(text, end, value);
    return parsed.ec == std::errc() && parsed.ptr == end;
}

void print_parallel_usage(const char* program) {
    std::cerr << "v4mer 2026-07-16 - parallel Klein V4 k-mer counter\n\n";
    std::cerr << "Usage: " << program
              << " <input> <k> <output.txt> [-t|--threads N]\n";
}

}  // namespace parallel_v4mer

int main(int argc, char* argv[]) {
    using namespace parallel_v4mer;
    if (argc < 4 || ((argc - 4) % 2) != 0) {
        print_parallel_usage(argv[0]);
        return 1;
    }

    size_t k = 0;
    if (!parse_size(argv[2], k) || k < 1 || k > 127) {
        std::cerr << "Error: k must be between 1 and 127\n";
        return 1;
    }

    size_t hardware_threads = std::thread::hardware_concurrency();
    if (hardware_threads == 0) hardware_threads = 1;
    size_t thread_count = std::min<size_t>(hardware_threads, 8);
    size_t memory_budget_mb = 0;
    for (int arg = 4; arg < argc; arg += 2) {
        const std::string option = argv[arg];
        size_t value = 0;
        if (arg + 1 >= argc || !parse_size(argv[arg + 1], value)) {
            std::cerr << "Error: option value must be an integer\n";
            return 1;
        }
        if (option == "-t" || option == "--threads") {
            if (value < 1 || value > 256) {
                std::cerr << "Error: thread count must be between 1 and 256\n";
                return 1;
            }
            thread_count = value;
        } else if (option == "--memory-budget-mb") {
            memory_budget_mb = value;
        } else {
            std::cerr << "Error: unknown option " << option << '\n';
            return 1;
        }
    }

    const std::string input_file = argv[1];
    const std::string output_file = argv[3];
    struct stat file_stat {};
    size_t file_size = 0;
    if (stat(input_file.c_str(), &file_stat) == 0 && file_stat.st_size > 0) {
        file_size = static_cast<size_t>(file_stat.st_size);
    }

    std::cerr << "v4mer 2026-07-16 - parallel Klein V4 k-mer counter\n";
    std::cerr << "===================================================\n";
    std::cerr << "Input: " << input_file << '\n';
    std::cerr << "K: " << k << '\n';
    std::cerr << "Output: " << output_file << '\n';

    try {
        const size_t words = (2 * k + 63) / 64;
        if (words == 1) {
            return run_parallel<1>(input_file, k, output_file, file_size,
                                   thread_count, memory_budget_mb);
        }
        if (words == 2) {
            return run_parallel<2>(input_file, k, output_file, file_size,
                                   thread_count, memory_budget_mb);
        }
        if (words == 3) {
            return run_parallel<3>(input_file, k, output_file, file_size,
                                   thread_count, memory_budget_mb);
        }
        return run_parallel<4>(input_file, k, output_file, file_size,
                               thread_count, memory_budget_mb);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}

