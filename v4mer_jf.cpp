/**
 * @file v4mer_jf.cpp
 * @brief Convert V4mer tabular output to a native Jellyfish 2.3.1 database.
 *
 * The writer is intentionally independent of the Jellyfish headers and
 * libraries.  It emits the binary/sorted layout consumed by Jellyfish 2.3.1.
 * The file stores the ordinary Jellyfish canonical projection; the four V4
 * channels remain available only in the original V4mer representation.
 */

#include <algorithm>
#include <array>
#include <charconv>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

constexpr unsigned MAX_K = 127;
constexpr unsigned WORDS = 4;
constexpr unsigned COUNTER_BYTES = 8;

using Words = std::array<uint64_t, WORDS>;

struct WordLess {
    bool operator()(const Words& left, const Words& right) const {
        for (size_t i = WORDS; i-- > 0;) {
            if (left[i] != right[i]) return left[i] < right[i];
        }
        return false;
    }
};

struct Record {
    Words key{};
    uint64_t count = 0;
    uint64_t position = 0;
};

[[noreturn]] void fail(const std::string& message) {
    throw std::runtime_error("v4mer-jf: " + message);
}

unsigned parse_k(const std::string& value) {
    unsigned k = 0;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), k);
    if (result.ec != std::errc() || result.ptr != value.data() + value.size() ||
        k == 0 || k > MAX_K) {
        fail("k must be in the range 1..127");
    }
    return k;
}

uint64_t parse_count(const std::string& value) {
    uint64_t count = 0;
    const auto result = std::from_chars(value.data(), value.data() + value.size(), count);
    if (result.ec != std::errc() || result.ptr != value.data() + value.size()) {
        fail("invalid count '" + value + "'");
    }
    return count;
}

int base_code(char base) {
    switch (base) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

std::string canonical_kmer(std::string_view input) {
    std::string reverse_complement;
    reverse_complement.resize(input.size());
    for (size_t i = 0; i < input.size(); ++i) {
        const int code = base_code(input[input.size() - 1 - i]);
        if (code < 0) fail("invalid DNA symbol in k-mer '" + std::string(input) + "'");
        reverse_complement[i] = "TGCA"[static_cast<size_t>(code)];
    }

    std::string result(input);
    for (char& base : result) {
        if (base_code(base) < 0) fail("invalid DNA symbol in k-mer '" + result + "'");
        base = static_cast<char>(std::toupper(static_cast<unsigned char>(base)));
    }
    return reverse_complement < result ? reverse_complement : result;
}

Words encode_kmer(std::string_view kmer) {
    Words words{};
    const unsigned total_bits = static_cast<unsigned>(kmer.size() * 2);
    for (size_t i = 0; i < kmer.size(); ++i) {
        const int code = base_code(kmer[i]);
        if (code < 0) fail("invalid DNA symbol in k-mer '" + std::string(kmer) + "'");
        const unsigned bit = total_bits - 2U - static_cast<unsigned>(2 * i);
        words[bit / 64U] |= static_cast<uint64_t>(code) << (bit % 64U);
    }
    return words;
}

uint64_t next_database_size(size_t records) {
    const size_t wanted = std::max<size_t>(2, records > (std::numeric_limits<size_t>::max() / 2)
                                               ? records
                                               : records * 2);
    uint64_t size = 2;
    while (size < wanted) {
        if (size > (std::numeric_limits<uint64_t>::max() >> 1))
            fail("too many records for a Jellyfish database");
        size <<= 1;
    }
    return size;
}

std::string json_header(unsigned k, uint64_t database_size) {
    const unsigned key_bits = k * 2;
    const unsigned rows = std::min<unsigned>(64, key_bits);
    return "{\"alignment\":8,\"canonical\":true,\"counter_len\":8,\"format\":\"binary/sorted\","
           "\"key_len\":" + std::to_string(key_bits) + ",\"matrix1\":{\"c\":" +
           std::to_string(key_bits) + ",\"identity\":true,\"r\":" + std::to_string(rows) +
           "},\"max_reprobe\":0,\"reprobes\":[0],\"size\":" +
           std::to_string(database_size) + "}";
}

void write_header(std::ofstream& output, const std::string& json) {
    constexpr size_t HEADER_DIGITS = 9;
    constexpr size_t ALIGNMENT = 8;
    size_t header_length = json.size();
    const size_t remainder = (HEADER_DIGITS + header_length) % ALIGNMENT;
    if (remainder != 0) header_length += ALIGNMENT - remainder;
    if (header_length > 999999999U) fail("Jellyfish header is too large");

    output << std::setw(static_cast<int>(HEADER_DIGITS)) << std::setfill('0')
           << header_length;
    output.write(json.data(), static_cast<std::streamsize>(json.size()));
    const std::string padding(header_length - json.size(), '\0');
    output.write(padding.data(), static_cast<std::streamsize>(padding.size()));
}

void write_key(std::ofstream& output, const Words& key, unsigned key_bytes) {
    output.write(reinterpret_cast<const char*>(key.data()), static_cast<std::streamsize>(key_bytes));
}

void write_counter(std::ofstream& output, uint64_t count) {
    // Jellyfish stores integer fields in the host byte order.  The supported
    // release target is the little-endian Jellyfish 2.3.1 build used on Linux.
    for (unsigned byte = 0; byte < COUNTER_BYTES; ++byte) {
        const char value = static_cast<char>((count >> (byte * 8U)) & 0xffU);
        output.write(&value, 1);
    }
}

void convert(const std::string& input_path, const std::string& output_path,
             unsigned requested_k) {
    const uint16_t endian_probe = 1;
    if (*reinterpret_cast<const unsigned char*>(&endian_probe) != 1)
        fail("only little-endian systems are supported");

    std::ifstream input(input_path);
    if (!input) fail("cannot open input '" + input_path + "'");

    std::map<Words, uint64_t, WordLess> counts;
    std::string line;
    unsigned k = requested_k;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;

        const size_t separator = line.find_first_of("\t ");
        if (separator == std::string::npos)
            fail("line " + std::to_string(line_number) + " must contain k-mer and count");
        const std::string kmer = line.substr(0, separator);
        size_t count_start = separator;
        while (count_start < line.size() &&
               (line[count_start] == '\t' || line[count_start] == ' ')) ++count_start;
        if (count_start == line.size())
            fail("line " + std::to_string(line_number) + " has no count");
        const size_t count_end = line.find_first_of("\t ", count_start);
        if (count_end != std::string::npos)
            fail("line " + std::to_string(line_number) + " has extra fields");

        if (k == 0) k = parse_k(std::to_string(kmer.size()));
        if (kmer.size() != k)
            fail("line " + std::to_string(line_number) + " has a k-mer of the wrong length");

        const uint64_t count = parse_count(line.substr(count_start));
        const Words encoded = encode_kmer(canonical_kmer(kmer));
        auto [entry, inserted] = counts.emplace(encoded, count);
        if (!inserted) {
            if (std::numeric_limits<uint64_t>::max() - entry->second < count)
                fail("count overflow while aggregating duplicate k-mers");
            entry->second += count;
        }
    }

    if (k == 0) fail("cannot infer k from an empty input; use --k");
    const unsigned key_bits = k * 2;
    const unsigned key_bytes = (key_bits + 7U) / 8U;
    const uint64_t database_size = next_database_size(counts.size());

    std::vector<Record> records;
    records.reserve(counts.size());
    const uint64_t mask = database_size - 1;
    for (const auto& [key, count] : counts)
        records.push_back(Record{key, count, key[0] & mask});
    std::sort(records.begin(), records.end(), [](const Record& left, const Record& right) {
        if (left.position != right.position) return left.position < right.position;
        return WordLess{}(left.key, right.key);
    });

    std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
    if (!output) fail("cannot open output '" + output_path + "'");
    write_header(output, json_header(k, database_size));
    for (const Record& record : records) {
        write_key(output, record.key, key_bytes);
        write_counter(output, record.count);
    }
    if (!output) fail("error writing output '" + output_path + "'");
}

void usage() {
    std::cerr << "Usage: v4mer-jf INPUT.tsv OUTPUT.jf [--k K]\n";
}

} // namespace

int main(int argc, char** argv) {
    if (argc < 3 || argc > 5) {
        usage();
        return 2;
    }
    try {
        unsigned requested_k = 0;
        if (argc == 4 || argc == 5) {
            if (std::string(argv[3]) != "--k" || argc != 5) {
                usage();
                return 2;
            }
            requested_k = parse_k(argv[4]);
        }
        convert(argv[1], argv[2], requested_k);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << "\n";
        return 1;
    }
}
