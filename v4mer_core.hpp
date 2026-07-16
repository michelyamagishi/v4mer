#pragma once

/**
 * @file v4mer_core.hpp
 * @brief Self-contained v4mer core with bulk V₄ count updates
 * @author Michel Eduardo Beleza Yamagishi
 * @version 1.1
 * @date 2026
 *
 * Uses Klein four-group V₄ = {I, R, C, RC} internally for efficient canonicalization,
 * then outputs in Jellyfish-compatible format (min of forward and reverse-complement).
 *
 * Supports FASTA/FASTQ formats (plain or gzip-compressed).
 *
 * Memory-optimized version with:
 * - Compact 8-bit counts (with an overflow table for high-count k-mers)
 * - Growable Robin Hood hash table with an 85% maximum load
 * - Lexicographic V4 canonicalization with transform-specific counts
 * - Buffered input and output
 *
 * Compile: g++ -std=c++17 -O3 -march=native -flto -o v4mer v4mer.cpp -lz
 */

#include <iostream>
#include <string>
#include <cstdint>
#include <algorithm>
#include <charconv>
#include <cerrno>
#include <cstdio>
#include <cstring>
#include <limits>
#include <sys/mman.h>
#include <sys/stat.h>
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
 * - Four 8-bit inline counts with a 64-bit overflow table for counts > 254
 * - Robin Hood hashing with automatic growth at 85% load
 * - Zero-filled anonymous mappings avoid an O(capacity) initialization pass
 *
 * Naturally aligned entry sizes:
 *   WORDS=1: 16 bytes (k ≤ 32)
 *   WORDS=2: 24 bytes (k ≤ 64)
 *   WORDS=3: 32 bytes (k ≤ 96)
 *   WORDS=4: 40 bytes (k ≤ 127)
 */
template<size_t WORDS>
class CompactHashTable {
public:
    static constexpr uint8_t OVERFLOW_MARKER = 255;
    static constexpr size_t MAX_LOAD_PERCENT = 85;

    // Keep 64-bit keys naturally aligned. Packing this structure to 1-byte
    // alignment saves four bytes for WORDS=1 but makes key accesses undefined
    // behavior (and faults on strict-alignment architectures).
    struct __attribute__((packed)) Entry {
        uint64_t kmer[WORDS];      // Klein canonical k-mer
        uint8_t count_I;           // Count for Identity transform (max 254, 255 = overflow)
        uint8_t count_R;           // Count for Reverse transform
        uint8_t count_C;           // Count for Complement transform
        uint8_t count_RC;          // Count for Reverse-Complement transform

        bool is_empty() const {
            return (count_I | count_R | count_C | count_RC) == 0;
        }
        void mark_empty() {
            count_I = count_R = count_C = count_RC = 0;
        }
    };

    // Overflow entry stores full 64-bit counts for channels exceeding 254
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
    size_t resize_threshold_;
    size_t load_percent_;
    size_t entry_count_;
    bool use_mmap_;

    // Overflow table for k-mers with counts > 254
    std::unordered_map<OverflowKey, OverflowCounts, OverflowKeyHash> overflow_table_;

public:
    explicit CompactHashTable(size_t expected_entries, size_t capacity_hint = 0,
                              size_t load_percent = MAX_LOAD_PERCENT)
        : entries_(nullptr), capacity_(0), resize_threshold_(0),
          load_percent_(load_percent), entry_count_(0), use_mmap_(false) {
        if (load_percent_ < 50 || load_percent_ > 95) {
            throw std::invalid_argument("Hash table load factor must be 50..95 percent");
        }

        size_t required_capacity = expected_entries;
        if (expected_entries > 0) {
            if (expected_entries > (std::numeric_limits<size_t>::max() - 99) /
                                   100) {
                throw std::length_error("Estimated k-mer table is too large");
            }
            required_capacity = (expected_entries * 100 + 99) /
                                load_percent_;
        }
        if (capacity_hint != 0) {
            capacity_ = std::max(capacity_hint, required_capacity);
            if (capacity_ < 16) capacity_ = 16;
        } else {
            // Without a global capacity hint, use power-of-two growth. A
            // supplied hint allows the caller to control aggregate capacity.
            capacity_ = 1;
            while (capacity_ < required_capacity) {
                if (capacity_ > std::numeric_limits<size_t>::max() / 2) {
                    throw std::length_error("Estimated k-mer table is too large");
                }
                capacity_ <<= 1;
            }
            if (capacity_ < 1024) capacity_ = 1024;
        }
        resize_threshold_ = capacity_ * load_percent_ / 100;
        if (resize_threshold_ == 0) resize_threshold_ = 1;
        entries_ = allocate_entries(capacity_, use_mmap_);
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
    inline static size_t compute_hash(const uint64_t* kmer_data) {
        uint64_t hash = kmer_data[0];
        for (size_t i = 1; i < WORDS; ++i) {
            uint64_t word = kmer_data[i] +
                0x9e3779b97f4a7c15ULL * static_cast<uint64_t>(i);
            word ^= word >> 33;
            word *= 0xff51afd7ed558ccdULL;
            word ^= word >> 33;
            hash ^= word;
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
        return slot >= hash ? slot - hash : capacity_ - hash + slot;
    }

    inline size_t bucket_index(size_t hash) const {
        return hash % capacity_;
    }

    inline size_t next_slot(size_t slot) const {
        return slot + 1 == capacity_ ? 0 : slot + 1;
    }

    /**
     * @brief Insert or increment k-mer count using Robin Hood hashing
     *
     * Robin Hood hashing maintains O(1) average lookup at high load
     * by "stealing" slots from entries with shorter probe distances.
     */
    void insert_or_increment(const uint64_t* kmer_data, Transform transform) {
        insert_impl<false>(kmer_data, transform, 1, compute_hash(kmer_data));
    }

    /**
     * @brief Insert a k-mer or add many observations to one V4 channel.
     *
     * The parallel counter coalesces repeated updates before taking a shard
     * lock. Applying the aggregate in one probe is essential for repetitive
     * sequences and preserves the compact/overflow count representation.
     */
    void insert_or_add(const uint64_t* kmer_data, Transform transform,
                       uint64_t amount) {
        if (amount == 0) return;
        if (amount == 1) {
            insert_impl<false>(kmer_data, transform, 1, compute_hash(kmer_data));
        } else {
            insert_impl<true>(kmer_data, transform, amount, compute_hash(kmer_data));
        }
    }

    void insert_or_add_hashed(const uint64_t* kmer_data, Transform transform,
                              uint64_t amount, size_t hash) {
        if (amount == 0) return;
        if (amount == 1) {
            insert_impl<false>(kmer_data, transform, 1, hash);
        } else {
            insert_impl<true>(kmer_data, transform, amount, hash);
        }
    }

private:
    template<bool BULK>
    void insert_impl(const uint64_t* kmer_data, Transform transform,
                     uint64_t amount, size_t hash) {
        size_t pos = bucket_index(hash);
        size_t dist = 0;

        // Entry to potentially insert (only populated if we need to insert new)
        Entry to_insert;
        bool inserting_new = false;
        __builtin_prefetch(&entries_[pos], 1, 3);

        while (true) {
            Entry& slot = entries_[pos];

            if (dist < 8) {
                __builtin_prefetch(&entries_[next_slot(next_slot(pos))], 1, 2);
            }

            if (slot.is_empty()) {
                if (inserting_new) {
                    // Place the displaced entry
                    slot = to_insert;
                } else {
                    // Insert new k-mer
                    copy_kmer(slot.kmer, kmer_data);
                    slot.count_I = slot.count_R = slot.count_C = slot.count_RC = 0;
                    if constexpr (BULK) {
                        add_compact_count(slot, transform, amount);
                    } else {
                        increment_compact_count(slot, transform);
                    }
                    ++entry_count_;
                }
                if (entry_count_ > resize_threshold_) grow();
                return;
            }

            if (!inserting_new && kmers_equal(slot.kmer, kmer_data)) {
                if constexpr (BULK) {
                    add_compact_count(slot, transform, amount);
                } else {
                    increment_compact_count(slot, transform);
                }
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
                    if constexpr (BULK) {
                        add_compact_count(to_insert, transform, amount);
                    } else {
                        increment_compact_count(to_insert, transform);
                    }
                    inserting_new = true;
                    ++entry_count_;
                }

                // Swap with current slot
                Entry tmp = slot;
                slot = to_insert;
                to_insert = tmp;
                dist = slot_dist;
            }

            pos = next_slot(pos);
            ++dist;

            // Safety check (should never trigger with proper load factor)
            if (dist > capacity_) {
                throw std::overflow_error("Hash table full (probe distance exceeded capacity)");

            }
        }
    }

    static Entry* allocate_entries(size_t capacity, bool& used_mmap) {
        if (capacity > std::numeric_limits<size_t>::max() / sizeof(Entry)) {
            throw std::length_error("K-mer table allocation is too large");
        }

        const size_t alloc_size = capacity * sizeof(Entry);
        void* allocation = mmap(nullptr, alloc_size, PROT_READ | PROT_WRITE,
                                MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
        if (allocation != MAP_FAILED) {
            used_mmap = true;
#ifdef MADV_HUGEPAGE
            (void)madvise(allocation, alloc_size, MADV_HUGEPAGE);
#endif
            // Anonymous mappings are already zero-filled. A zero count tuple is
            // the empty marker, so no O(capacity) initialization pass is needed.
            return static_cast<Entry*>(allocation);
        }

        used_mmap = false;
        Entry* entries = new Entry[capacity];
        std::memset(entries, 0, alloc_size);
        return entries;
    }

    static void release_entries(Entry* entries, size_t capacity, bool used_mmap) {
        if (used_mmap) {
            (void)munmap(entries, capacity * sizeof(Entry));
        } else {
            delete[] entries;
        }
    }

    void insert_existing(Entry entry) {
        size_t hash = compute_hash(entry.kmer);
        size_t pos = bucket_index(hash);
        size_t dist = 0;

        while (true) {
            Entry& slot = entries_[pos];
            if (slot.is_empty()) {
                slot = entry;
                return;
            }

            size_t slot_hash = compute_hash(slot.kmer);
            size_t slot_dist = probe_distance(pos, slot_hash);
            if (slot_dist < dist) {
                std::swap(slot, entry);
                dist = slot_dist;
            }
            pos = next_slot(pos);
            ++dist;
        }
    }

    void grow() {
        if (capacity_ > std::numeric_limits<size_t>::max() / 2) {
            throw std::length_error("K-mer table cannot grow further");
        }

        Entry* old_entries = entries_;
        const size_t old_capacity = capacity_;
        const bool old_use_mmap = use_mmap_;
        const size_t new_capacity = capacity_ * 2;
        bool new_use_mmap = false;
        Entry* new_entries = allocate_entries(new_capacity, new_use_mmap);

        capacity_ = new_capacity;
        resize_threshold_ = capacity_ * load_percent_ / 100;
        entries_ = new_entries;
        use_mmap_ = new_use_mmap;

        for (size_t i = 0; i < old_capacity; ++i) {
            if (!old_entries[i].is_empty()) insert_existing(old_entries[i]);
        }
        release_entries(old_entries, old_capacity, old_use_mmap);
    }

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
            ++(*count_ptr);
        } else {
            add_overflow(entry.kmer, transform, 1);
        }
    }

    inline void add_compact_count(Entry& entry, Transform transform,
                                  uint64_t amount) {
        uint8_t* count_ptr = nullptr;
        switch (transform) {
            case TRANSFORM_I:  count_ptr = &entry.count_I;  break;
            case TRANSFORM_R:  count_ptr = &entry.count_R;  break;
            case TRANSFORM_C:  count_ptr = &entry.count_C;  break;
            case TRANSFORM_RC: count_ptr = &entry.count_RC; break;
        }

        if (*count_ptr < OVERFLOW_MARKER) {
            const uint64_t room = OVERFLOW_MARKER - *count_ptr;
            if (amount <= room) {
                *count_ptr = static_cast<uint8_t>(*count_ptr + amount);
                return;
            }
            *count_ptr = OVERFLOW_MARKER;
            amount -= room;
        }
        add_overflow(entry.kmer, transform, amount);
    }

    void add_overflow(const uint64_t* kmer_data, Transform transform,
                      uint64_t amount) {
        OverflowKey key;
        copy_kmer(key.kmer, kmer_data);

        auto& counts = overflow_table_[key];
        uint64_t* target = nullptr;
        switch (transform) {
            case TRANSFORM_I:  target = &counts.count_I;  break;
            case TRANSFORM_R:  target = &counts.count_R;  break;
            case TRANSFORM_C:  target = &counts.count_C;  break;
            case TRANSFORM_RC: target = &counts.count_RC; break;
        }
        if (amount > std::numeric_limits<uint64_t>::max() - *target) {
            throw std::overflow_error("K-mer count exceeds uint64_t");
        }
        *target += amount;
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
    size_t capacity() const { return capacity_; }
    double load_factor() const {
        return static_cast<double>(entry_count_) /
               static_cast<double>(capacity_);
    }

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
            const uint64_t differing_bits = data_[i] ^ other.data_[i];
            if (differing_bits != 0) {
                const unsigned bit =
                    static_cast<unsigned>(__builtin_ctzll(differing_bits)) & ~1U;
                const uint8_t lhs = static_cast<uint8_t>((data_[i] >> bit) & 3);
                const uint8_t rhs =
                    static_cast<uint8_t>((other.data_[i] >> bit) & 3);
                return lhs < rhs;
            }
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
        for (size_t i = 0; i < WORDS; ++i) {
            const uint64_t differing_bits = a[i] ^ b[i];
            if (differing_bits != 0) {
                const unsigned bit =
                    static_cast<unsigned>(__builtin_ctzll(differing_bits)) & ~1U;
                const uint8_t lhs = static_cast<uint8_t>((a[i] >> bit) & 3);
                const uint8_t rhs = static_cast<uint8_t>((b[i] >> bit) & 3);
                return lhs < rhs ? -1 : 1;
            }
        }
        return 0;
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

        // Complement reverses the A<C<G<T alphabet. Therefore F versus C(F)
        // is decided by the first base alone, as is RC(F) versus R(F). This
        // leaves one full lexicographic comparison instead of three.
        CompactKmer<WORDS> complement_of_F;
        const uint64_t* left;
        Transform left_transform;
        if ((forward_.data()[0] & 3) < 2) {
            left = forward_.data();
            left_transform = TRANSFORM_I;
        } else {
            complement_of_F = compute_complement(forward_);
            left = complement_of_F.data();
            left_transform = TRANSFORM_C;
        }

        CompactKmer<WORDS> reverse_of_F;
        const uint64_t* right;
        Transform right_transform;
        if ((reverse_complement_.data()[0] & 3) < 2) {
            right = reverse_complement_.data();
            right_transform = TRANSFORM_RC;
        } else {
            reverse_of_F = compute_complement(reverse_complement_);
            right = reverse_of_F.data();
            right_transform = TRANSFORM_R;
        }

        const auto precedence = [](Transform transform) -> unsigned {
            switch (transform) {
                case TRANSFORM_I:  return 0;
                case TRANSFORM_RC: return 1;
                case TRANSFORM_C:  return 2;
                case TRANSFORM_R:  return 3;
            }
            return 4;
        };

        const int comparison = compare_kmers(right, left);
        if (comparison < 0 ||
            (comparison == 0 && precedence(right_transform) <
                                    precedence(left_transform))) {
            result.transform = right_transform;
        } else {
            result.transform = left_transform;
        }
        const uint64_t* best =
            (result.transform == right_transform) ? right : left;

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
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;
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
                int error_number = Z_OK;
                const char* message = gzerror(gz_file_, &error_number);
                throw std::runtime_error(std::string("Error reading gzip file: ") +
                                         (message ? message : "unknown error"));
            }
            valid_ = static_cast<size_t>(bytes_read);
            if (bytes_read == 0) {
                eof_ = true;
            }
        } else {
            size_t bytes_read = fread(buffer_.data(), 1, BUFFER_SIZE, plain_file_);
            valid_ = bytes_read;
            if (bytes_read == 0) {
                if (ferror(plain_file_)) {
                    throw std::runtime_error("Error reading input file");
                }
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
    FILE* file = fopen(filename.c_str(), "rb");
    if (!file) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    unsigned char magic[2] = {0, 0};
    const size_t bytes_read = fread(magic, 1, sizeof(magic), file);
    const bool read_error = ferror(file) != 0;
    fclose(file);
    if (read_error) throw std::runtime_error("Cannot read file: " + filename);

    return bytes_read == sizeof(magic) && magic[0] == 0x1f && magic[1] == 0x8b
               ? CompressionType::GZIP
               : CompressionType::NONE;
}

FileFormat detect_format(BufferedReader& reader) {
    int c = reader.peek();
    while (c == '\n' || c == '\r' || c == ' ' || c == '\t') {
        reader.get();
        c = reader.peek();
    }
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
                if (c == '>') window.reset();
                continue;
            }

            at_line_start = false;
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

    template<size_t WORDS>
    void count_fastq(BufferedReader& reader, size_t k, CompactHashTable<WORDS>& table) {
        CompactRollingWindow<WORDS> window(k);

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
            window.reset();
            size_t sequence_length = 0;

            // FASTQ permits wrapped sequence and quality lines. Continue sequence
            // lines until a '+' separator appears at the beginning of a line.
            while (true) {
                c = reader.peek();
                if (c == EOF) {
                    throw std::runtime_error("Malformed FASTQ: missing '+' line");
                }
                if (c == '+') {
                    reader.skip_line();
                    break;
                }

                while ((c = reader.get()) != '\n' && c != EOF) {
                    if (c == '\r') continue;
                    if (sequence_length == std::numeric_limits<size_t>::max()) {
                        throw std::length_error("FASTQ sequence is too long");
                    }
                    ++sequence_length;
                    int8_t code = BASE_ENCODING[static_cast<unsigned char>(c)];
                    if (code >= 0) {
                        if (window.add_base(static_cast<uint8_t>(code))) {
                            auto canonical = window.get_canonical();
                            table.insert_or_increment(canonical.kmer.data(),
                                                      canonical.transform);
                        }
                    } else {
                        window.reset();
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
                    throw std::runtime_error("Malformed FASTQ: truncated quality data");
                }
                if (c != '\n' && c != '\r') ++quality_length;
            }

            // The final quality line must end exactly at sequence_length.
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
    }
};

// ============================================================================
// UTILITY FUNCTIONS
// ============================================================================

inline uint8_t encoded_base(const uint64_t* kmer_data, size_t position) {
    const size_t bit_position = 2 * position;
    return static_cast<uint8_t>(
        (kmer_data[bit_position >> 6] >> (bit_position & 63)) & 3);
}

enum class OutputTransform { IDENTITY, REVERSE, COMPLEMENT };

class BufferedOutput {
private:
    static constexpr size_t BUFFER_SIZE = 1024 * 1024;
    FILE* file_ = nullptr;
    std::string buffer_;

    void flush() {
        if (buffer_.empty()) return;
        const size_t bytes_written =
            fwrite(buffer_.data(), 1, buffer_.size(), file_);
        if (bytes_written != buffer_.size()) {
            throw std::runtime_error("Error writing output file");
        }
        buffer_.clear();
    }

public:
    explicit BufferedOutput(const std::string& filename) {
        file_ = fopen(filename.c_str(), "wb");
        if (!file_) {
            throw std::runtime_error("Cannot open output file: " + filename +
                                     ": " + std::strerror(errno));
        }
        buffer_.reserve(BUFFER_SIZE + 256);
    }

    ~BufferedOutput() {
        if (file_) fclose(file_);
    }

    BufferedOutput(const BufferedOutput&) = delete;
    BufferedOutput& operator=(const BufferedOutput&) = delete;

    void write(const uint64_t* kmer_data, size_t k, OutputTransform transform,
               uint64_t count) {
        for (size_t i = 0; i < k; ++i) {
            const size_t source = transform == OutputTransform::REVERSE
                                      ? k - 1 - i
                                      : i;
            uint8_t base = encoded_base(kmer_data, source);
            if (transform == OutputTransform::COMPLEMENT) base ^= 3;
            buffer_.push_back(BASE_DECODING[base]);
        }
        buffer_.push_back('\t');

        char count_buffer[32];
        const auto conversion =
            std::to_chars(count_buffer, count_buffer + sizeof(count_buffer), count);
        if (conversion.ec != std::errc()) {
            throw std::runtime_error("Cannot format k-mer count");
        }
        buffer_.append(count_buffer, conversion.ptr);
        buffer_.push_back('\n');

        if (buffer_.size() >= BUFFER_SIZE) flush();
    }

    void close() {
        flush();
        FILE* file = file_;
        file_ = nullptr;
        if (fclose(file) != 0) {
            throw std::runtime_error("Error closing output file");
        }
    }
};

bool reverse_is_canonical_pair_member(const uint64_t* kmer_data, size_t k) {
    for (size_t i = 0; i < k; ++i) {
        const uint8_t reverse_base = encoded_base(kmer_data, k - 1 - i);
        const uint8_t complement_base = encoded_base(kmer_data, i) ^ 3;
        if (reverse_base != complement_base) {
            return reverse_base < complement_base;
        }
    }
    return true;
}

/**
 * @brief Estimate unique k-mers based on file characteristics
 *
 * Uses format-aware heuristics to avoid over-allocation:
 * - FASTQ sequence is at most about half of a valid file (quality has equal length)
 * - Gzip typically compresses genomic data 3-5x
 * - Caps estimate at 4^k (maximum possible unique k-mers)
 */
size_t estimate_unique_kmers(size_t file_size, size_t k, FileFormat format, CompressionType compression) {
    size_t seq_bytes = file_size;

    if (compression == CompressionType::GZIP) {
        seq_bytes = file_size > std::numeric_limits<size_t>::max() / 4
                        ? std::numeric_limits<size_t>::max()
                        : file_size * 4;
    }

    if (format == FileFormat::FASTQ) {
        // A valid FASTQ contains sequence and quality of the same length,
        // plus headers and separators. Reads also commonly repeat across
        // sequencing cycles, so planning at the raw-base count overallocates
        // the table. Keep a conservative 30% unique-class estimate; the hash
        // table still grows safely if an unusual input exceeds it.
        seq_bytes = (seq_bytes / 5) * 2;
        seq_bytes = (seq_bytes / 10) * 3;
    }

    const size_t max_possible =
        k <= 31 ? (static_cast<size_t>(1) << (2 * k))
                : std::numeric_limits<size_t>::max();
    return std::min(seq_bytes, max_possible);
}

#ifndef V4MER_CORE_ONLY

void print_usage(const char* program) {
    std::cerr << "v4mer 1.1 - Klein V₄ k-mer counter (Jellyfish-compatible)\n\n";
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
    std::cerr << "Using " << WORDS << "-word entries ("
              << sizeof(typename CompactHashTable<WORDS>::Entry)
              << " bytes each)\n\n";

    // Detect format first for better estimation
    SequenceReader reader(input_file);

    const char* format_str = (reader.format() == FileFormat::FASTA) ? "FASTA" : "FASTQ";
    const char* comp_str = (reader.compression() == CompressionType::GZIP) ? " (gzip)" : "";
    std::cerr << "Format: " << format_str << comp_str << "\n";

    // Use improved capacity estimation
    size_t estimated_kmers = estimate_unique_kmers(file_size, k, reader.format(), reader.compression());
    std::cerr << "Estimated unique k-mers: " << estimated_kmers << "\n";

    CompactHashTable<WORDS> table(estimated_kmers);

    std::cerr << "Counting canonical k-mers...\n";

    reader.count_kmers(k, table);

    std::cerr << "Found " << table.size() << " distinct equivalence classes\n";
    std::cerr << "Hash table capacity: " << table.capacity() << " entries\n";
    std::cerr << "Hash table load factor: " << (table.load_factor() * 100) << "%\n";
    if (table.overflow_size() > 0) {
        std::cerr << "Overflow entries (counts > 254): " << table.overflow_size() << "\n";
    }
    std::cerr << "Writing output...\n";

    BufferedOutput output(output_file);

    size_t output_lines = 0;

    for (const auto& entry : table) {
        CompactKmer<WORDS> canonical;
        for (size_t i = 0; i < WORDS; ++i) {
            canonical.data()[i] = entry.kmer[i];
        }

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
            // The V4 representative is the minimum of all four transforms, so
            // it is necessarily canonical for its {I, RC} Jellyfish pair.
            output.write(canonical.data(), k, OutputTransform::IDENTITY,
                         count_pair1);
            ++output_lines;
        }

        if (count_pair2 > 0) {
            const OutputTransform pair_transform =
                reverse_is_canonical_pair_member(canonical.data(), k)
                    ? OutputTransform::REVERSE
                    : OutputTransform::COMPLEMENT;
            output.write(canonical.data(), k, pair_transform, count_pair2);
            ++output_lines;
        }
    }
    output.close();

    std::cerr << "Wrote " << output_lines << " output lines\n";
    std::cerr << "Done.\n";
    return 0;
}

int main(int argc, char* argv[]) {
    if (argc != 4) {
        print_usage(argv[0]);
        return 1;
    }

    size_t k = 0;
    const char* k_begin = argv[2];
    const char* k_end = k_begin + std::strlen(k_begin);
    const auto parsed_k = std::from_chars(k_begin, k_end, k);
    if (parsed_k.ec != std::errc() || parsed_k.ptr != k_end || k < 1 || k > 127) {
        std::cerr << "Error: k must be between 1 and 127\n";
        return 1;
    }

    const std::string input_file = argv[1];
    const std::string output_file = argv[3];

    std::cerr << "v4mer 1.1 - Klein V₄ k-mer counter\n";
    std::cerr << "==================================\n";
    std::cerr << "Input:  " << input_file << "\n";
    std::cerr << "K:      " << k << "\n";
    std::cerr << "Output: " << output_file << "\n";

    struct stat file_stat;
    size_t file_size = 0;
    if (stat(input_file.c_str(), &file_stat) == 0) {
        file_size = static_cast<size_t>(file_stat.st_size);
    }

    try {
        const size_t words_needed = (2 * k + 63) / 64;
        if (words_needed == 1) {
            return run_counting<1>(input_file, k, output_file, file_size);
        }
        if (words_needed == 2) {
            return run_counting<2>(input_file, k, output_file, file_size);
        }
        if (words_needed == 3) {
            return run_counting<3>(input_file, k, output_file, file_size);
        }
        return run_counting<4>(input_file, k, output_file, file_size);
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}

#endif  // V4MER_CORE_ONLY
