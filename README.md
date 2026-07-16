# v4mer

> Parallel DNA k-mer counter based on full Klein four-group (V₄) canonicalization.

V4mer is designed for reproducible, memory-controlled k-mer counting with
Jellyfish-compatible output. It is distributed as a single executable and has
no runtime dependency beyond zlib.

## At a glance

| Capability | Description |
|---|---|
| Canonicalization | Full V₄: identity (I), reverse (R), complement (C), and reverse-complement (RC) |
| Input | FASTA, FASTQ, and gzip-compressed files |
| k range | `1..127` |
| Parallelism | `1..256` worker threads |
| Memory control | Optional `--memory-budget-mb` limit |
| Temporary files | None; counting remains in memory |
| Output | Tab-separated k-mer/count pairs compatible with Jellyfish |

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

### Options

| Option | Meaning | Default |
|---|---|---:|
| `-t N`, `--threads N` | Number of worker threads; valid range `1..256` | Up to 8 hardware threads |
| `--memory-budget-mb N` | Maximum planned memory in MiB; execution is rejected before counting if exceeded | `0` (automatic budget) |

## Input and output

Input format is detected automatically. Ambiguous bases reset the current k-mer
window. Wrapped FASTA, CRLF line endings, and wrapped FASTQ records are
supported.

Each output line contains a canonical k-mer and its count:

```text
AAAAAAAAAAAAAAAAAAAAA	1523
AAAAAAAAAAAAAAAAAAAAT	892
```

For direct comparison with Jellyfish, normalize both outputs to
`kmer<TAB>count`, sort by k-mer and count, and compare the resulting streams.
The validated release outputs contain the same k-mers with the same counts as
Jellyfish for every tested dataset and k value.

## Benchmark results

The benchmark used four real inputs, `k=21,31,33,57,101`, 1/2/4/8 threads,
and two repetitions. Timing includes counting and text output for both tools.
Jellyfish used canonical mode (`-C`) followed by text dumping.

### Overall results

| Metric | Result |
|---|---:|
| Complete dataset/k validations | 20 |
| Median benchmark cases faster than Jellyfish | 80 / 80 |
| Minimum end-to-end speedup | 1.46x |
| Maximum V4mer/Jellyfish RSS ratio | 1.96x |
| Minimum 4-thread counting speedup vs. 1 thread | 2.73x |
| Worst 8-thread end-to-end change vs. 4 threads | +0.8% |

### Representative end-to-end results at four threads

| Dataset | Format | k | V4mer (s) | Jellyfish (s) | Speedup | RSS ratio |
|---|---|---:|---:|---:|---:|---:|
| Arcopilus | FASTA | 101 | 6.350 | 23.945 | 3.77x | 0.90x |
| HuRef chr22 | FASTA | 101 | 6.280 | 17.700 | 2.82x | 0.90x |
| MT147 | FASTQ | 57 | 3.670 | 7.720 | 2.10x | 0.42x |
| MT147 | FASTQ gzip | 101 | 4.910 | 8.365 | 1.70x | 0.34x |

`RSS ratio` is V4mer peak resident memory divided by Jellyfish peak resident
memory for the same dataset, k, and thread count.

### Scaling examples

| Dataset | Format | k | 1 thread (s) | 4 threads (s) | 8 threads (s) | 8 vs. 4 |
|---|---|---:|---:|---:|---:|---:|
| Arcopilus | FASTA | 101 | 17.245 | 6.350 | 5.855 | -7.8% |
| HuRef chr22 | FASTA | 101 | 16.960 | 6.280 | 6.040 | -3.8% |
| MT147 | FASTQ | 57 | 11.510 | 3.670 | 3.240 | -11.7% |
| MT147 | FASTQ | 101 | 11.755 | 4.445 | 4.480 | +0.8% |
| MT147 | FASTQ gzip | 101 | 12.655 | 4.910 | 4.555 | -7.2% |

The complete measurements were performed before packaging. V4mer/Jellyfish
validation compares normalized key/count streams; output order and whitespace
are not used as cross-tool equality criteria. V4mer output is deterministic
across its own thread counts.

Results depend on CPU, compiler, storage, input composition, and Jellyfish
configuration. Gzip decompression and physical output writes remain serial.

## License

MIT License. See [LICENSE](LICENSE).

## Citation

Yamagishi, Michel Eduardo Beleza. *Mathematical Grammar of Biology*.
SpringerBriefs in Mathematics. Cham: Springer, 2017. DOI:
[10.1007/978-3-319-62689-5](https://doi.org/10.1007/978-3-319-62689-5).
eBook ISBN: 978-3-319-62689-5. Softcover ISBN: 978-3-319-62688-8.

```bibtex
@book{yamagishi2017mathematical,
  author    = {Yamagishi, Michel Eduardo Beleza},
  title     = {Mathematical Grammar of Biology},
  series    = {SpringerBriefs in Mathematics},
  publisher = {Springer},
  address   = {Cham},
  year      = {2017},
  doi       = {10.1007/978-3-319-62689-5},
  isbn      = {978-3-319-62689-5}
}
```
