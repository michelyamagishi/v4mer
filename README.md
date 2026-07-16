# v4mer

> Parallel DNA k-mer counter based on full Klein four-group (V₄) canonicalization.

V4mer is designed for reproducible, memory-controlled k-mer counting with
Jellyfish-compatible output. It is distributed as two small executables and has
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
| Output | Deterministic TSV; optional native Jellyfish `.jf` conversion |

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

The resulting executables are `./v4mer` and `./v4mer-jf`.

## Usage

```bash
./v4mer INPUT K OUTPUT [OPTIONS]
```

Examples:

```bash
./v4mer genome.fa 21 kmers.txt --threads 8
./v4mer reads.fastq.gz 31 kmers.txt --threads 8 --memory-budget-mb 4096
./v4mer-jf kmers.txt kmers.jf
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

### Native Jellyfish output

Some downstream tools require Jellyfish's binary database format rather than
text. Convert the deterministic V4mer output with the separate `v4mer-jf`
utility:

```bash
./v4mer INPUT K OUTPUT.tsv --threads 8
./v4mer-jf OUTPUT.tsv OUTPUT.jf
jellyfish info OUTPUT.jf
jellyfish query OUTPUT.jf ACGT...
jellyfish dump -c -t OUTPUT.jf
```

`v4mer-jf` is self-contained and targets the native `binary/sorted` format of
Jellyfish 2.3.1 on little-endian systems. The converter canonicalizes and
aggregates the same ordinary `{I,RC}` key/count projection consumed by
Jellyfish. A `.jf` file does not preserve the four independent V₄ channels;
retain the V4mer TSV output when those channels are required.

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

### Complete four-thread comparison

All values below are median end-to-end times from the measured matrix.

| Dataset | Format | k | V4mer (s) | Jellyfish (s) | Speedup | RSS ratio |
|---|---|---:|---:|---:|---:|---:|
| Arcopilus | FASTA | 21 | 2.300 | 6.455 | 2.81x | 1.95x |
| Arcopilus | FASTA | 31 | 2.640 | 8.060 | 3.05x | 1.16x |
| Arcopilus | FASTA | 33 | 2.855 | 7.330 | 2.57x | 1.81x |
| Arcopilus | FASTA | 57 | 3.870 | 9.750 | 2.52x | 0.99x |
| Arcopilus | FASTA | 101 | 6.350 | 23.945 | 3.77x | 0.90x |
| HuRef chr22 | FASTA | 21 | 2.020 | 5.530 | 2.74x | 1.93x |
| HuRef chr22 | FASTA | 31 | 2.485 | 6.715 | 2.70x | 1.17x |
| HuRef chr22 | FASTA | 33 | 2.830 | 7.075 | 2.50x | 1.82x |
| HuRef chr22 | FASTA | 57 | 3.710 | 9.540 | 2.57x | 0.97x |
| HuRef chr22 | FASTA | 101 | 6.280 | 17.700 | 2.82x | 0.90x |
| MT147 | FASTQ | 21 | 2.710 | 6.265 | 2.31x | 0.85x |
| MT147 | FASTQ | 31 | 2.900 | 6.785 | 2.34x | 0.51x |
| MT147 | FASTQ | 33 | 3.165 | 6.955 | 2.20x | 0.78x |
| MT147 | FASTQ | 57 | 3.670 | 7.720 | 2.10x | 0.42x |
| MT147 | FASTQ | 101 | 4.445 | 9.020 | 2.03x | 0.38x |
| MT147 | FASTQ gzip | 21 | 3.150 | 6.215 | 1.97x | 0.75x |
| MT147 | FASTQ gzip | 31 | 3.400 | 6.690 | 1.97x | 0.46x |
| MT147 | FASTQ gzip | 33 | 3.705 | 7.460 | 2.01x | 0.69x |
| MT147 | FASTQ gzip | 57 | 4.450 | 7.985 | 1.79x | 0.43x |
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
