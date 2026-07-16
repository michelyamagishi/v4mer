# v4mer

A parallel DNA k-mer counter based on full Klein four-group (V₄) canonicalization.
This release is the parallel implementation published as `v4mer`.

## Theoretical framework

V4mer is based on the V₄ symmetry framework described in *Mathematical Grammar
of Biology*. For fixed-length DNA words, the identity (I), reverse (R),
complement (C), and reverse-complement (RC) transformations form the Klein
four-group, V₄ = {I, R, C, RC}. These transformations partition k-mers into
symmetry-equivalence classes.

V4mer selects a canonical representative for each class and retains the four
transform-specific count channels. This extends the usual forward/reverse-
complement pairing while preserving Jellyfish-compatible k-mer/count output.

Reference: Michel Eduardo Beleza Yamagishi, *Mathematical Grammar of Biology*,
Springer, 2017. [Springer book page](https://link.springer.com/book/10.1007/978-3-319-62689-5).

## Features

- Full V₄ canonicalization: identity, reverse, complement, and reverse-complement.
- Jellyfish-compatible tab-separated k-mer/count output.
- Deterministic output across thread counts.
- Supports `k=1..127`.
- Supports FASTA, FASTQ, and gzip-compressed input.
- Bounded in-memory execution with a configurable memory budget.
- No temporary counting files or spill-to-disk.

## Build

Requirements: a C++17 compiler, GNU Make, and zlib.

```bash
make
```

The resulting executable is `./v4mer`.

## Usage

```bash
./v4mer INPUT K OUTPUT [OPTIONS]
```

Examples:

```bash
./v4mer genome.fa 21 kmers.txt --threads 8
./v4mer reads.fastq.gz 31 kmers.txt --threads 8 --memory-budget-mb 4096
```

Options:

| Option | Description |
|---|---|
| `-t N`, `--threads N` | Use `N` worker threads; valid range: 1--256. |
| `--memory-budget-mb N` | Reject execution before counting if the planned memory exceeds `N` MiB. `0` selects the automatic budget. |

Input format is detected automatically. Ambiguous bases reset the current k-mer
window. FASTA wrapping, CRLF line endings, and wrapped FASTQ records are supported.

## Output and canonicalization

Each window is assigned to the lexicographically smallest member of its V₄
equivalence class. The four transform channels remain separate internally and
are emitted in Jellyfish-compatible key/count form.

For a direct comparison with Jellyfish, normalize each output to
`kmer<TAB>count`, sort by k-mer and count, and compare the resulting streams.
The measured release outputs contain the same k-mers with the same counts as
Jellyfish for every validated dataset and k value.

## Benchmark results

The benchmark used four real inputs, `k=21,31,33,57,101`, 1/2/4/8 threads,
and two repetitions. Timing includes counting and text output for both tools.
Jellyfish used canonical mode (`-C`) followed by text dumping.

| Metric | Result |
|---|---:|
| Complete dataset/k validations | 20 |
| V4mer faster than Jellyfish | 80/80 median cases |
| Minimum end-to-end speedup | 1.46x |
| Maximum V4mer/Jellyfish RSS ratio | 1.96x |
| Minimum 4-thread counting speedup vs. 1 thread | 2.73x |
| Worst 8-thread end-to-end change vs. 4 threads | +0.8% |

Representative median results at four threads:

| Input | k | V4mer (s) | Jellyfish (s) | Speedup |
|---|---:|---:|---:|---:|
| Arcopilus FASTA | 101 | 6.350 | 23.945 | 3.77x |
| HuRef FASTA | 101 | 6.280 | 17.700 | 2.82x |
| FASTQ plain | 57 | 3.670 | 7.720 | 2.10x |
| FASTQ gzip | 101 | 4.910 | 8.365 | 1.70x |

Results depend on CPU, compiler, storage, input composition, and Jellyfish
configuration. Gzip decompression and physical output writes remain serial.

## License

MIT License. See [LICENSE](LICENSE).
