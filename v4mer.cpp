/**
 * @file v4mer.cpp
 * @brief v4mer - Klein V₄ k-mer counter (Jellyfish-compatible output)
 * @author Michel Eduardo Beleza Yamagishi
 * @version 1.0
 * @date 2026
 *
 * Uses Klein four-group V₄ = {I, R, C, RC} internally for efficient canonicalization,
 * then outputs in Jellyfish-compatible format (min of forward and reverse-complement).
 *
 * Supports FASTA/FASTQ formats (plain or gzip-compressed).
 *
 * Memory-optimized version with:
 * - Compact 16-bit counts (with overflow table for rare high-count k-mers)
 * - Robin Hood hashing for 90% load factor
 * - Better capacity estimation
 * - Zero-copy canonical k-mer computation
 *
 * Compile: g++ -std=c++17 -O3 -march=native -flto -o v4mer v4mer.cpp -lz
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <algorithm>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <zlib.h>
#include <vector>
#include <stdexcept>
#include <unordered_map>

// ============================================================================
// FILE FORMAT AND COMPRESSION ENUMS
// ============================================================================

enum class FileFormat { FASTA, FASTQ };
enum class CompressionType { NONE, GZIP };

// ============================================================================
// CONSTANTS AND LOOKUP TABLES
// ============================================================================

static constexpr int8_t BASE_ENCODING[256] = {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  0, -1,  1, -1, -1, -1,  2, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
};

static constexpr char BASE_DECODING[4] = {'A', 'C', 'G', 'T'};

// Transform indices for Klein four-group
enum Transform : uint8_t {
    TRANSFORM_I  = 0,  // Identity (forward)
    TRANSFORM_R  = 1,  // Reverse
    TRANSFORM_C  = 2,  // Complement
    TRANSFORM_RC = 3   // Reverse-Complement
};

// ============================================================================
// TEMPLATE-BASED COMPACT HASH TABLE WITH ROBIN HOOD HASHING
// ============================================================================

/**
 * @brief Memory-optimized hash table with Robin Hood hashing
 * @tparam WORDS Number of 64-bit words needed for k-mer storage
 *
 * Optimizations:
 * - 8-bit counts (4 bytes vs 8 bytes for 4 counts) - 25% smaller entries
 * - Robin Hood hashing allows 90% load factor (vs 50% with linear probing)
 * - Overflow table for rare k-mers with counts > 254
 *
 * Entry sizes (with packing):
 *   WORDS=1: 12 bytes (k ≤ 31)  - was 16 bytes (25% reduction)
 *   WORDS=2: 20 bytes (k ≤ 63)  - was 24 bytes (17% reduction)
 *   WORDS=3: 28 bytes (k ≤ 95)  - was 32 bytes (12% reduction)
 *   WORDS=4: 36 bytes (k ≤ 127) - was 40 bytes (10% reduction)
 */
template<size_t WORDS>
class CompactHashTable {
public:
    static constexpr uint64_t EMPTY_MARKER = ~0ULL;
    static constexpr uint8_t OVERFLOW_MARKER = 255;

    // Pack the struct to avoid padding between kmer and counts
    #pragma pack(push, 1)
    struct Entry {
        uint64_t kmer[WORDS];      // Klein canonical k-mer
        uint8_t count_I;           // Count for Identity transform (max 254, 255 = overflow)
        uint8_t count_R;           // Count for Reverse transform
        uint8_t count_C;           // Count for Complement transform
        uint8_t count_RC;          // Count for Reverse-Complement transform

        bool is_empty() const { return kmer[0] == EMPTY_MARKER; }
        void mark_empty() {
            kmer[0] = EMPTY_MARKER;
            count_I = count_R = count_C = count_RC = 0;
        }
    };
    #pragma pack(pop)

    // Overflow entry stores full 64-bit counts for k-mers exceeding 16-bit limit
    struct OverflowKey {
        uint64_t kmer[WORDS];

        bool operator==(const OverflowKey& other) const {
            for (size_t i = 0; i < WORDS; ++i) {
                if (kmer[i] != other.kmer[i]) return false;
            }
            return true;
        }
    };

    struct OverflowKeyHash {
        size_t operator()(const OverflowKey& key) const {
            size_t hash = key.kmer[0];
            for (size_t i = 1; i < WORDS; ++i) {
                hash ^= key.kmer[i];
            }
            hash ^= hash >> 33;
            hash *= 0xff51afd7ed558ccdULL;
            hash ^= hash >> 33;
            return hash;
        }
    };

    struct OverflowCounts {
        uint64_t count_I = 0;
        uint64_t count_R = 0;
        uint64_t count_C = 0;
        uint64_t count_RC = 0;
    };

private:
    Entry* entries_;
    size_t capacity_;
    size_t mask_;
    size_t entry_count_;
    size_t kmer_length_;
    bool use_mmap_;

    // Overflow table for k-mers with counts > 254
    std::unordered_map<OverflowKey, OverflowCounts, OverflowKeyHash> overflow_table_;

public:
    CompactHashTable(size_t k, size_t expected_entries)
        : entry_count_(0), kmer_length_(k), use_mmap_(false) {

        // Size table as power of 2, targeting ~90% load factor (Robin Hood allows this)
        capacity_ = 1;
        while (capacity_ < static_cast<size_t>(expected_entries * 1.15)) {
            capacity_ <<= 1;
        }
        if (capacity_ < 1024) capacity_ = 1024;
        mask_ = capacity_ - 1;

        // Use mmap for large allocations (better memory management)
        size_t alloc_size = capacity_ * sizeof(Entry);
        if (alloc_size > 1024 * 1024) {
            entries_ = static_cast<Entry*>(mmap(nullptr, alloc_size,
                PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0));
            if (entries_ != MAP_FAILED) {
                use_mmap_ = true;
                madvise(entries_, alloc_size, MADV_HUGEPAGE);
            } else {
                entries_ = new Entry[capacity_];
            }
        } else {
            entries_ = new Entry[capacity_];
        }

        // Initialize all entries as empty
        for (size_t i = 0; i < capacity_; ++i) {
            entries_[i].mark_empty();
        }
    }

    ~CompactHashTable() {
        if (use_mmap_) {
            munmap(entries_, capacity_ * sizeof(Entry));
        } else {
            delete[] entries_;
        }
    }

    CompactHashTable(const CompactHashTable&) = delete;
    CompactHashTable& operator=(const CompactHashTable&) = delete;

    /**
     * @brief Fast hash function (Murmur3 64-bit finalizer)
     */
    inline size_t compute_hash(const uint64_t* kmer_data) const {
        size_t hash = kmer_data[0];
        for (size_t i = 1; i < WORDS; ++i) {
            hash ^= kmer_data[i];
        }
        hash ^= hash >> 33;
        hash *= 0xff51afd7ed558ccdULL;
        hash ^= hash >> 33;
        hash *= 0xc4ceb9fe1a85ec53ULL;
        hash ^= hash >> 33;
        return hash;
    }

    /**
     * @brief Compare k-mers for equality (unrolled for small WORDS)
     */
    inline bool kmers_equal(const uint64_t* a, const uint64_t* b) const {
        if constexpr (WORDS == 1) {
            return a[0] == b[0];
        } else if constexpr (WORDS == 2) {
            return a[0] == b[0] && a[1] == b[1];
        } else {
            for (size_t i = 0; i < WORDS; ++i) {
                if (a[i] != b[i]) return false;
            }
            return true;
        }
    }

    /**
     * @brief Compute probe distance for Robin Hood hashing
     */
    inline size_t probe_distance(size_t slot, size_t hash) const {
        return (slot - hash) & mask_;
    }

    /**
     * @brief Insert or increment k-mer count using Robin Hood hashing
     *
     * Robin Hood hashing maintains O(1) average lookup even at 90%+ load
     * by "stealing" slots from entries with shorter probe distances.
     */
    void insert_or_increment(const uint64_t* kmer_data, Transform transform) {
        size_t hash = compute_hash(kmer_data);
        size_t pos = hash & mask_;
        size_t dist = 0;

        // Entry to potentially insert (only populated if we need to insert new)
        Entry to_insert;
        bool inserting_new = false;
        Transform insert_transform = transform;

        __builtin_prefetch(&entries_[pos], 1, 3);

        while (true) {
            Entry& slot = entries_[pos];

            if (dist < 8) {
                __builtin_prefetch(&entries_[(pos + 2) & mask_], 1, 2);
            }

            if (slot.is_empty()) {
                if (inserting_new) {
                    // Place the displaced entry
                    slot = to_insert;
                } else {
                    // Insert new k-mer
                    copy_kmer(slot.kmer, kmer_data);
                    slot.count_I = slot.count_R = slot.count_C = slot.count_RC = 0;
                    increment_compact_count(slot, transform);
                    ++entry_count_;
                }
                return;
            }

            if (!inserting_new && kmers_equal(slot.kmer, kmer_data)) {
                // Found existing k-mer, increment count
                increment_compact_count(slot, transform);
                return;
            }

            // Robin Hood: check if we should swap
            size_t slot_hash = compute_hash(slot.kmer);
            size_t slot_dist = probe_distance(pos, slot_hash);

            if (slot_dist < dist) {
                // Current entry has traveled less than us - swap (Robin Hood)
                if (!inserting_new) {
                    // First time swapping - create entry for the new k-mer
                    copy_kmer(to_insert.kmer, kmer_data);
                    to_insert.count_I = to_insert.count_R = to_insert.count_C = to_insert.count_RC = 0;
                    increment_compact_count(to_insert, transform);
                    inserting_new = true;
                    ++entry_count_;
                }

                // Swap with current slot
                Entry tmp = slot;
                slot = to_insert;
                to_insert = tmp;
                dist = slot_dist;
            }

            pos = (pos + 1) & mask_;
            ++dist;

            // Safety check (should never trigger with proper load factor)
            if (dist > capacity_) {
                std::cerr << "Error: Hash table full (probe distance exceeded capacity)\n";
                exit(1);
            }
        }
    }

private:
    inline void copy_kmer(uint64_t* dst, const uint64_t* src) {
        if constexpr (WORDS == 1) {
            dst[0] = src[0];
        } else if constexpr (WORDS == 2) {
            dst[0] = src[0];
            dst[1] = src[1];
        } else {
            for (size_t i = 0; i < WORDS; ++i) {
                dst[i] = src[i];
            }
        }
    }

    inline void increment_compact_count(Entry& entry, Transform transform) {
        uint8_t* count_ptr = nullptr;
        switch (transform) {
            case TRANSFORM_I:  count_ptr = &entry.count_I;  break;
            case TRANSFORM_R:  count_ptr = &entry.count_R;  break;
            case TRANSFORM_C:  count_ptr = &entry.count_C;  break;
            case TRANSFORM_RC: count_ptr = &entry.count_RC; break;
        }

        if (*count_ptr < OVERFLOW_MARKER) {
            // Normal case: increment 8-bit counter (0 -> 254)
            ++(*count_ptr);
        } else {
            // Already at overflow marker (255) - use overflow table
            increment_overflow(entry.kmer, transform);
        }
    }

    void increment_overflow(const uint64_t* kmer_data, Transform transform) {
        OverflowKey key;
        copy_kmer(key.kmer, kmer_data);

        auto& counts = overflow_table_[key];
        switch (transform) {
            case TRANSFORM_I:  ++counts.count_I;  break;
            case TRANSFORM_R:  ++counts.count_R;  break;
            case TRANSFORM_C:  ++counts.count_C;  break;
            case TRANSFORM_RC: ++counts.count_RC; break;
        }
    }

public:
    /**
     * @brief Get the full count for a transform, handling overflow
     */
    uint64_t get_count(const Entry& entry, Transform transform) const {
        uint8_t compact;
        switch (transform) {
            case TRANSFORM_I:  compact = entry.count_I;  break;
            case TRANSFORM_R:  compact = entry.count_R;  break;
            case TRANSFORM_C:  compact = entry.count_C;  break;
            case TRANSFORM_RC: compact = entry.count_RC; break;
            default: return 0;
        }

        if (compact < OVERFLOW_MARKER) {
            return compact;
        }

        // Look up in overflow table (count >= 255)
        OverflowKey key;
        for (size_t i = 0; i < WORDS; ++i) {
            key.kmer[i] = entry.kmer[i];
        }

        auto it = overflow_table_.find(key);
        if (it == overflow_table_.end()) {
            return OVERFLOW_MARKER;  // Shouldn't happen, but safe fallback
        }

        // Return 255 (marker) + overflow count
        switch (transform) {
            case TRANSFORM_I:  return OVERFLOW_MARKER + it->second.count_I;
            case TRANSFORM_R:  return OVERFLOW_MARKER + it->second.count_R;
            case TRANSFORM_C:  return OVERFLOW_MARKER + it->second.count_C;
            case TRANSFORM_RC: return OVERFLOW_MARKER + it->second.count_RC;
            default: return 0;
        }
    }

    size_t size() const { return entry_count_; }
    size_t overflow_size() const { return overflow_table_.size(); }
    double load_factor() const { return static_cast<double>(entry_count_) / capacity_; }

    // Iterator support
    class Iterator {
        const CompactHashTable* table_;
        size_t index_;
        void skip_empty() {
            while (index_ < table_->capacity_ && table_->entries_[index_].is_empty()) ++index_;
        }
    public:
        Iterator(const CompactHashTable* t, size_t i) : table_(t), index_(i) { skip_empty(); }
        bool operator!=(const Iterator& o) const { return index_ != o.index_; }
        Iterator& operator++() { ++index_; skip_empty(); return *this; }
        const Entry& operator*() const { return table_->entries_[index_]; }
    };

    Iterator begin() const { return Iterator(this, 0); }
    Iterator end() const { return Iterator(this, capacity_); }
};

// ============================================================================
// COMPACT PACKED KMER (Template-based)
// ============================================================================

template<size_t WORDS>
class CompactKmer {
private:
    uint64_t data_[WORDS];

public:
    CompactKmer() { clear(); }

    void clear() {
        if constexpr (WORDS == 1) {
            data_[0] = 0;
        } else if constexpr (WORDS == 2) {
            data_[0] = 0;
            data_[1] = 0;
        } else {
            for (size_t i = 0; i < WORDS; ++i) data_[i] = 0;
        }
    }

    uint64_t* data() { return data_; }
    const uint64_t* data() const { return data_; }

    bool operator<(const CompactKmer& other) const {
        for (size_t i = 0; i < WORDS; ++i) {
            if (data_[i] != other.data_[i]) return data_[i] < other.data_[i];
        }
        return false;
    }
};

/**
 * @brief Result of canonical form computation
 * Contains both the canonical k-mer and which transform was applied
 */
template<size_t WORDS>
struct CanonicalResult {
    CompactKmer<WORDS> kmer;
    Transform transform;
};

// ============================================================================
// COMPACT DUAL ROLLING WINDOW
// ============================================================================

template<size_t WORDS>
class CompactRollingWindow {
private:
    size_t k_;
    CompactKmer<WORDS> forward_;
    CompactKmer<WORDS> reverse_complement_;
    size_t bases_in_window_;
    uint64_t high_word_mask_;
    uint64_t complement_mask_[WORDS];  // Mask for XOR complement

public:
    explicit CompactRollingWindow(size_t k) : k_(k), bases_in_window_(0) {
        size_t used_bits = (2 * k) % 64;
        high_word_mask_ = (used_bits == 0) ? ~0ULL : ((1ULL << used_bits) - 1);

        // Initialize complement mask: all 1s in the 2*k bits used for k-mer storage
        // XOR with this mask computes DNA complement: A(00)↔T(11), C(01)↔G(10)
        size_t total_bits = 2 * k;
        for (size_t w = 0; w < WORDS; ++w) {
            size_t bits_in_word = (total_bits >= 64) ? 64 : total_bits;
            complement_mask_[w] = (bits_in_word == 64) ? ~0ULL : ((1ULL << bits_in_word) - 1);
            total_bits = (total_bits > 64) ? total_bits - 64 : 0;
        }
    }

    void reset() {
        forward_.clear();
        reverse_complement_.clear();
        bases_in_window_ = 0;
    }

    bool add_base(uint8_t base) {
        uint8_t complement = 3 - base;

        if (bases_in_window_ < k_) {
            size_t fwd_bit = 2 * bases_in_window_;
            size_t fwd_word = fwd_bit >> 6;
            size_t fwd_offset = fwd_bit & 63;
            forward_.data()[fwd_word] |= static_cast<uint64_t>(base) << fwd_offset;

            size_t rc_bit = 2 * (k_ - 1 - bases_in_window_);
            size_t rc_word = rc_bit >> 6;
            size_t rc_offset = rc_bit & 63;
            reverse_complement_.data()[rc_word] |= static_cast<uint64_t>(complement) << rc_offset;

            ++bases_in_window_;
            return bases_in_window_ == k_;
        }

        uint64_t* fwd = forward_.data();
        uint64_t* rc = reverse_complement_.data();

        // Forward: shift right by 2
        if constexpr (WORDS == 1) {
            fwd[0] >>= 2;
        } else {
            for (size_t i = 0; i < WORDS; ++i) {
                fwd[i] >>= 2;
                if (i + 1 < WORDS) fwd[i] |= (fwd[i + 1] & 3) << 62;
            }
        }

        size_t high_bit = 2 * (k_ - 1);
        size_t high_word = high_bit >> 6;
        size_t high_offset = high_bit & 63;
        fwd[high_word] = (fwd[high_word] & ~(3ULL << high_offset)) | (static_cast<uint64_t>(base) << high_offset);
        fwd[WORDS - 1] &= high_word_mask_;

        // Reverse complement: shift left by 2
        if constexpr (WORDS == 1) {
            rc[0] <<= 2;
            rc[0] = (rc[0] & ~3ULL) | complement;
        } else {
            for (size_t i = WORDS; i > 0; --i) {
                size_t idx = i - 1;
                rc[idx] <<= 2;
                if (idx > 0) rc[idx] |= (rc[idx - 1] >> 62) & 3;
            }
            rc[0] = (rc[0] & ~3ULL) | complement;
        }
        rc[WORDS - 1] &= high_word_mask_;

        return true;
    }

    /**
     * @brief Compute DNA complement using XOR
     *
     * DNA complement swaps bases: A(00) ↔ T(11), C(01) ↔ G(10)
     * This is equivalent to XOR with 0b11 for each 2-bit base,
     * i.e., XOR with a mask of all 1s in the used bit positions.
     *
     * Time complexity: O(WORDS) - very fast!
     */
    CompactKmer<WORDS> compute_complement(const CompactKmer<WORDS>& kmer) const {
        CompactKmer<WORDS> result;
        if constexpr (WORDS == 1) {
            result.data()[0] = kmer.data()[0] ^ complement_mask_[0];
        } else if constexpr (WORDS == 2) {
            result.data()[0] = kmer.data()[0] ^ complement_mask_[0];
            result.data()[1] = kmer.data()[1] ^ complement_mask_[1];
        } else {
            for (size_t i = 0; i < WORDS; ++i) {
                result.data()[i] = kmer.data()[i] ^ complement_mask_[i];
            }
        }
        return result;
    }

    /**
     * @brief Compare two k-mers lexicographically
     * @return negative if a < b, 0 if equal, positive if a > b
     */
    inline int compare_kmers(const uint64_t* a, const uint64_t* b) const {
        if constexpr (WORDS == 1) {
            return (a[0] < b[0]) ? -1 : ((a[0] > b[0]) ? 1 : 0);
        } else if constexpr (WORDS == 2) {
            if (a[0] != b[0]) return (a[0] < b[0]) ? -1 : 1;
            return (a[1] < b[1]) ? -1 : ((a[1] > b[1]) ? 1 : 0);
        } else {
            for (size_t i = 0; i < WORDS; ++i) {
                if (a[i] != b[i]) return (a[i] < b[i]) ? -1 : 1;
            }
            return 0;
        }
    }

    /**
     * @brief Get canonical form and transform using Klein four-group V₄ (zero-copy optimized)
     *
     * Uses the Klein four-group V₄ = {I, R, C, RC} to find the canonical form.
     * Optimized to minimize copies: computes C and R lazily, only copies the winner.
     *
     * @return CanonicalResult with canonical k-mer and applied transform
     */
    CanonicalResult<WORDS> get_canonical() const {
        CanonicalResult<WORDS> result;

        // Start with forward (I) as candidate
        const uint64_t* best = forward_.data();
        result.transform = TRANSFORM_I;

        // Compare against RC (already maintained)
        if (compare_kmers(reverse_complement_.data(), best) < 0) {
            best = reverse_complement_.data();
            result.transform = TRANSFORM_RC;
        }

        // Compute C = complement(F) and compare
        CompactKmer<WORDS> complement_of_F = compute_complement(forward_);
        if (compare_kmers(complement_of_F.data(), best) < 0) {
            best = complement_of_F.data();
            result.transform = TRANSFORM_C;
        }

        // Compute R = complement(RC) and compare
        CompactKmer<WORDS> reverse_of_F = compute_complement(reverse_complement_);
        if (compare_kmers(reverse_of_F.data(), best) < 0) {
            best = reverse_of_F.data();
            result.transform = TRANSFORM_R;
        }

        // Copy only the winner
        if constexpr (WORDS == 1) {
            result.kmer.data()[0] = best[0];
        } else if constexpr (WORDS == 2) {
            result.kmer.data()[0] = best[0];
            result.kmer.data()[1] = best[1];
        } else {
            for (size_t i = 0; i < WORDS; ++i) {
                result.kmer.data()[i] = best[i];
            }
        }

        return result;
    }
};

// ============================================================================
// BUFFERED READER (Supports plain and gzip files)
// ============================================================================

class BufferedReader {
private:
    static constexpr size_t BUFFER_SIZE = 16 * 1024 * 1024;  // 16MB
    std::vector<char> buffer_;
    size_t pos_ = 0;
    size_t valid_ = 0;
    gzFile gz_file_ = nullptr;
    FILE* plain_file_ = nullptr;
    CompressionType compression_;
    bool eof_ = false;

    void fill_buffer() {
        if (eof_) return;

        if (compression_ == CompressionType::GZIP) {
            int bytes_read = gzread(gz_file_, buffer_.data(), BUFFER_SIZE);
            if (bytes_read < 0) {
                throw std::runtime_error("Error reading gzip file");
            }
            valid_ = static_cast<size_t>(bytes_read);
            if (bytes_read == 0) {
                eof_ = true;
            }
        } else {
            size_t bytes_read = fread(buffer_.data(), 1, BUFFER_SIZE, plain_file_);
            valid_ = bytes_read;
            if (bytes_read == 0) {
                eof_ = true;
            }
        }
        pos_ = 0;
    }

public:
    BufferedReader(const std::string& filename, CompressionType comp)
        : buffer_(BUFFER_SIZE), compression_(comp) {

        if (compression_ == CompressionType::GZIP) {
            gz_file_ = gzopen(filename.c_str(), "rb");
            if (!gz_file_) {
                throw std::runtime_error("Cannot open gzip file: " + filename);
            }
            gzbuffer(gz_file_, BUFFER_SIZE);
        } else {
            plain_file_ = fopen(filename.c_str(), "rb");
            if (!plain_file_) {
                throw std::runtime_error("Cannot open file: " + filename);
            }
        }

        fill_buffer();
    }

    ~BufferedReader() {
        if (gz_file_) gzclose(gz_file_);
        if (plain_file_) fclose(plain_file_);
    }

    BufferedReader(const BufferedReader&) = delete;
    BufferedReader& operator=(const BufferedReader&) = delete;

    int peek() {
        if (pos_ >= valid_) {
            fill_buffer();
            if (eof_) return EOF;
        }
        return static_cast<unsigned char>(buffer_[pos_]);
    }

    int get() {
        if (pos_ >= valid_) {
            fill_buffer();
            if (eof_) return EOF;
        }
        return static_cast<unsigned char>(buffer_[pos_++]);
    }

    void skip_line() {
        int c;
        while ((c = get()) != '\n' && c != EOF) {
            // Skip until newline or EOF
        }
    }

    bool eof() const { return eof_ && pos_ >= valid_; }
};

// ============================================================================
// FORMAT/COMPRESSION DETECTION
// ============================================================================

CompressionType detect_compression(const std::string& filename) {
    if (filename.size() >= 3 && filename.substr(filename.size() - 3) == ".gz") {
        return CompressionType::GZIP;
    }
    return CompressionType::NONE;
}

FileFormat detect_format(BufferedReader& reader) {
    int c = reader.peek();
    if (c == '>') return FileFormat::FASTA;
    if (c == '@') return FileFormat::FASTQ;
    throw std::runtime_error("Unknown file format: expected '>' (FASTA) or '@' (FASTQ)");
}

// ============================================================================
// SEQUENCE READER (Unified FASTA/FASTQ interface)
// ============================================================================

class SequenceReader {
private:
    std::string filename_;
    CompressionType compression_;
    FileFormat format_;

public:
    explicit SequenceReader(const std::string& filename)
        : filename_(filename) {
        compression_ = detect_compression(filename);

        // Open briefly to detect format
        BufferedReader reader(filename, compression_);
        format_ = detect_format(reader);
    }

    FileFormat format() const { return format_; }
    CompressionType compression() const { return compression_; }

    template<size_t WORDS>
    void count_kmers(size_t k, CompactHashTable<WORDS>& table) {
        BufferedReader reader(filename_, compression_);

        if (format_ == FileFormat::FASTA) {
            count_fasta<WORDS>(reader, k, table);
        } else {
            count_fastq<WORDS>(reader, k, table);
        }
    }

private:
    template<size_t WORDS>
    void count_fasta(BufferedReader& reader, size_t k, CompactHashTable<WORDS>& table) {
        CompactRollingWindow<WORDS> window(k);
        bool in_header = false;

        int c;
        while ((c = reader.get()) != EOF) {
            if (c == '>') {
                in_header = true;
                window.reset();
            } else if (c == '\n') {
                in_header = false;
            } else if (!in_header) {
                int8_t code = BASE_ENCODING[static_cast<unsigned char>(c)];
                if (code >= 0) {
                    if (window.add_base(static_cast<uint8_t>(code))) {
                        auto canonical = window.get_canonical();
                        table.insert_or_increment(canonical.kmer.data(), canonical.transform);
                    }
                } else {
                    window.reset();
                }
            }
        }
    }

    template<size_t WORDS>
    void count_fastq(BufferedReader& reader, size_t k, CompactHashTable<WORDS>& table) {
        CompactRollingWindow<WORDS> window(k);

        // FASTQ format: 4 lines per record
        // Line 1: @header
        // Line 2: sequence
        // Line 3: +
        // Line 4: quality scores

        while (!reader.eof()) {
            int c = reader.peek();
            if (c == EOF) break;

            // Skip any blank lines or unexpected characters
            if (c != '@') {
                reader.skip_line();
                continue;
            }

            // Line 1: Skip header line
            reader.skip_line();

            // Line 2: Process sequence line
            while ((c = reader.get()) != '\n' && c != EOF) {
                int8_t code = BASE_ENCODING[static_cast<unsigned char>(c)];
                if (code >= 0) {
                    if (window.add_base(static_cast<uint8_t>(code))) {
                        auto canonical = window.get_canonical();
                        table.insert_or_increment(canonical.kmer.data(), canonical.transform);
                    }
                } else {
                    window.reset();
                }
            }

            // Line 3: Skip + line
            reader.skip_line();

            // Line 4: Skip quality line
            reader.skip_line();

            // Reset window between reads (each read is independent)
            window.reset();
        }
    }
};

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

std::string kmer_to_string(const uint64_t* kmer_data, size_t k) {
    std::string result;
    result.reserve(k);
    for (size_t i = 0; i < k; ++i) {
        size_t bit_pos = 2 * i;
        size_t word = bit_pos >> 6;
        size_t offset = bit_pos & 63;
        uint8_t base = (kmer_data[word] >> offset) & 3;
        result += BASE_DECODING[base];
    }
    return result;
}

/**
 * @brief Compute complement of a k-mer (XOR with all-1s mask)
 */
template<size_t WORDS>
CompactKmer<WORDS> compute_complement_kmer(const CompactKmer<WORDS>& kmer, size_t k) {
    CompactKmer<WORDS> result;
    size_t total_bits = 2 * k;
    for (size_t w = 0; w < WORDS; ++w) {
        size_t bits_in_word = (total_bits >= 64) ? 64 : total_bits;
        uint64_t mask = (bits_in_word == 64) ? ~0ULL : ((1ULL << bits_in_word) - 1);
        result.data()[w] = kmer.data()[w] ^ mask;
        total_bits = (total_bits > 64) ? total_bits - 64 : 0;
    }
    return result;
}

/**
 * @brief Compute reverse of a k-mer (reverse the order of bases)
 * O(k) operation, only used during output
 */
template<size_t WORDS>
CompactKmer<WORDS> compute_reverse_kmer(const CompactKmer<WORDS>& kmer, size_t k) {
    CompactKmer<WORDS> result;
    for (size_t i = 0; i < k; ++i) {
        // Get base at position i
        size_t src_bit = 2 * i;
        size_t src_word = src_bit >> 6;
        size_t src_offset = src_bit & 63;
        uint8_t base = (kmer.data()[src_word] >> src_offset) & 3;

        // Put at position (k-1-i)
        size_t dst_bit = 2 * (k - 1 - i);
        size_t dst_word = dst_bit >> 6;
        size_t dst_offset = dst_bit & 63;
        result.data()[dst_word] |= static_cast<uint64_t>(base) << dst_offset;
    }
    return result;
}

/**
 * @brief Compute reverse-complement of a k-mer
 */
template<size_t WORDS>
CompactKmer<WORDS> compute_revcomp_kmer(const CompactKmer<WORDS>& kmer, size_t k) {
    return compute_complement_kmer(compute_reverse_kmer(kmer, k), k);
}

/**
 * @brief Estimate unique k-mers based on file characteristics
 *
 * Uses format-aware heuristics to avoid over-allocation:
 * - FASTQ has ~4x overhead (header, +, quality lines)
 * - Gzip typically compresses genomic data 3-5x
 * - Caps estimate at 4^k (maximum possible unique k-mers)
 */
size_t estimate_unique_kmers(size_t file_size, size_t k, FileFormat format, CompressionType compression) {
    // Estimate actual sequence bytes
    size_t seq_bytes = file_size;

    // For gzipped files, estimate decompressed size (typically 3-5x larger)
    if (compression == CompressionType::GZIP) {
        seq_bytes = file_size * 4;  // Conservative estimate
    }

    // FASTQ: ~25% is sequence (header, seq, +, quality = 4 lines, seq is ~1/4)
    // FASTA: ~90% is sequence (just headers to skip)
    if (format == FileFormat::FASTQ) {
        seq_bytes = seq_bytes / 4;
    } else {
        seq_bytes = (seq_bytes * 9) / 10;
    }

    // Maximum possible unique k-mers is 4^k, but cap at 2^62 to avoid overflow
    size_t max_possible = (k <= 31) ? (1ULL << (2 * k)) : (1ULL << 62);

    // For genomic data with k >= 15, most k-mers are unique
    // Total possible k-mers in sequence: seq_bytes - k + 1
    // Empirically, Klein V4 equivalence classes ≈ seq_bytes for bacterial genomes
    // Use 120% of seq_bytes as safety margin (observed: actual can exceed estimate)
    size_t estimated = (seq_bytes * 12) / 10;

    // Minimum reasonable size
    if (estimated < 1000000) estimated = 1000000;

    return std::min(estimated, max_possible);
}

void print_usage(const char* program) {
    std::cerr << "v4mer 1.0 - Klein V₄ k-mer counter (Jellyfish-compatible)\n\n";
    std::cerr << "Usage: " << program << " <input> <k> <output.txt>\n\n";
    std::cerr << "Supported formats (auto-detected):\n";
    std::cerr << "  FASTA:  .fa, .fasta, .fa.gz, .fasta.gz\n";
    std::cerr << "  FASTQ:  .fq, .fastq, .fq.gz, .fastq.gz\n\n";
    std::cerr << "Examples:\n";
    std::cerr << "  " << program << " genome.fa 29 output.txt\n";
    std::cerr << "  " << program << " reads.fastq.gz 21 output.txt\n";
}

// ============================================================================
// MAIN - Dispatches to correct template instantiation based on k
// ============================================================================

template<size_t WORDS>
int run_counting(const std::string& input_file, size_t k, const std::string& output_file,
                 size_t file_size) {
    std::cerr << "Using " << WORDS << "-word entries (" << (WORDS * 8 + 4) << " bytes each)\n\n";

    // Detect format first for better estimation
    SequenceReader reader(input_file);

    const char* format_str = (reader.format() == FileFormat::FASTA) ? "FASTA" : "FASTQ";
    const char* comp_str = (reader.compression() == CompressionType::GZIP) ? " (gzip)" : "";
    std::cerr << "Format: " << format_str << comp_str << "\n";

    // Use improved capacity estimation
    size_t estimated_kmers = estimate_unique_kmers(file_size, k, reader.format(), reader.compression());
    std::cerr << "Estimated unique k-mers: " << estimated_kmers << "\n";

    CompactHashTable<WORDS> table(k, estimated_kmers);

    std::cerr << "Counting canonical k-mers...\n";

    reader.count_kmers(k, table);

    std::cerr << "Found " << table.size() << " distinct equivalence classes\n";
    std::cerr << "Hash table load factor: " << (table.load_factor() * 100) << "%\n";
    if (table.overflow_size() > 0) {
        std::cerr << "Overflow entries (counts > 254): " << table.overflow_size() << "\n";
    }
    std::cerr << "Writing output...\n";

    std::ofstream output(output_file);
    if (!output) {
        std::cerr << "Error: Cannot open output file\n";
        return 1;
    }

    size_t output_lines = 0;

    for (const auto& entry : table) {
        CompactKmer<WORDS> canonical;
        for (size_t i = 0; i < WORDS; ++i) {
            canonical.data()[i] = entry.kmer[i];
        }

        // Compute all 4 variants from the Klein canonical
        CompactKmer<WORDS> kmer_C = compute_complement_kmer(canonical, k);
        CompactKmer<WORDS> kmer_R = compute_reverse_kmer(canonical, k);
        CompactKmer<WORDS> kmer_RC = compute_revcomp_kmer(canonical, k);

        // Get full counts (handles overflow automatically)
        uint64_t count_I = table.get_count(entry, TRANSFORM_I);
        uint64_t count_R = table.get_count(entry, TRANSFORM_R);
        uint64_t count_C = table.get_count(entry, TRANSFORM_C);
        uint64_t count_RC = table.get_count(entry, TRANSFORM_RC);

        // Jellyfish pair 1: {I, RC} - forward and reverse-complement
        // Jellyfish pair 2: {R, C}  - reverse and complement
        uint64_t count_pair1 = count_I + count_RC;
        uint64_t count_pair2 = count_R + count_C;

        if (count_pair1 > 0) {
            CompactKmer<WORDS> jf_canon1 = (canonical < kmer_RC) ? canonical : kmer_RC;
            output << kmer_to_string(jf_canon1.data(), k) << '\t' << count_pair1 << '\n';
            ++output_lines;
        }

        if (count_pair2 > 0) {
            CompactKmer<WORDS> jf_canon2 = (kmer_R < kmer_C) ? kmer_R : kmer_C;
            output << kmer_to_string(jf_canon2.data(), k) << '\t' << count_pair2 << '\n';
            ++output_lines;
        }
    }

    std::cerr << "Wrote " << output_lines << " output lines\n";
    std::cerr << "Done.\n";
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }

    std::string input_file = argv[1];
    int k = std::atoi(argv[2]);
    std::string output_file = argv[3];

    if (k < 1 || k > 127) {
        std::cerr << "Error: k must be between 1 and 127\n";
        return 1;
    }

    std::cerr << "v4mer 1.0 - Klein V₄ k-mer counter\n";
    std::cerr << "==================================\n";
    std::cerr << "Input:  " << input_file << "\n";
    std::cerr << "K:      " << k << "\n";
    std::cerr << "Output: " << output_file << "\n";

    struct stat file_stat;
    size_t file_size = 0;
    if (stat(input_file.c_str(), &file_stat) == 0) {
        file_size = static_cast<size_t>(file_stat.st_size);
    }

    size_t words_needed = (2 * k + 63) / 64;

    if (words_needed == 1) {
        return run_counting<1>(input_file, k, output_file, file_size);
    } else if (words_needed == 2) {
        return run_counting<2>(input_file, k, output_file, file_size);
    } else if (words_needed == 3) {
        return run_counting<3>(input_file, k, output_file, file_size);
    } else {
        return run_counting<4>(input_file, k, output_file, file_size);
    }
}
