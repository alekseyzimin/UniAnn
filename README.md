# UniAnn

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

UniAnn is a universal **ab initio gene finder** for annotating protein-coding genes in genome sequences. It combines nucleotide-state emissions from a coding-frame model with start, stop, donor, and acceptor site scores, then uses a Viterbi hidden Markov model to produce a GFF3 gene annotation.


## Features

- Predicts CDS, intron, and intergenic regions from one FASTA sequence
- Uses a seven-state Viterbi model:
  - `N`: noncoding/intergenic
  - `E0`, `E1`, `E2`: exon states in three reading frames
  - `I0`, `I1`, `I2`: intron states in three frames
- Enforces ATG starts, TAA/TAG/TGA stops, GT donor sites, and AG acceptor sites
- Enforces minimum exon, intron, intergenic, and single-exon lengths
- Produces GFF3 output through bundled `gffread`
- Includes a Perl preprocessing script for per-base coding-frame scores
- Includes C++ code for computing splice-site scores from zoeHMM-style positive and negative models
- Includes a source distribution and installation script

## Repository layout

```text
UniAnn/
├── src/
│   ├── Makefile
│   ├── uniann.cpp                 # Viterbi gene decoder
│   └── compute_markov_scores.cpp  # zoeHMM splice-site score calculator
├── scripts/
│   ├── uniann.sh                  # main annotation pipeline
│   ├── convert_coding_frame.pl    # coding-frame score converter
│   └── uniann_both_strands.pl     # optional experimental both-strand helper
├── gffcompare                     # bundled gffcompare executable
├── gffread                        # bundled gffread executable
├── install.sh                     # builds and installs UniAnn into ./bin
├── make_dist.sh                   # creates uniann.1.0.0.tgz
├── LICENSE                       # GNU GPL v3
└── uniann.1.0.0.tgz              # bundled source distribution
```

## Requirements

### System tools

- Linux or another Unix-like environment
- Bash
- GNU Make
- A C++17 compiler, such as `g++`
- Perl
- `gffread` and `gffcompare` (both are included in the repository)

The Makefile builds with `-std=c++17 -O3 -Wall -Wextra -march=native`.

### Input data

The main pipeline requires:

1. A FASTA file containing **one sequence**
2. A tab-delimited coding-frame score file
3. A tab-delimited training/annotation score file containing start, stop, donor, and acceptor probabilities

The repository does not bundle example input FASTA or model files. You must provide compatible inputs, or prepare them from the training model used with UniAnn.

## Installation

Clone the repository and run the installer from the repository root:

```bash
git clone https://github.com/colindaven/UniAnn.git
cd UniAnn
./install.sh
```

The installer:

1. Creates `./bin`
2. Copies `gffcompare` and `gffread` into `./bin`
3. Builds `uniann` and `compute_markov_scores` in `src/`
4. Moves the C++ executables to `./bin`
5. Copies Perl and Bash scripts to `./bin`

After installation, add UniAnn to your `PATH` if desired:

```bash
export PATH="$PWD/bin:$PATH"
```

The binaries and scripts in `bin/` are not stored in Git; they are generated locally by `install.sh`.

### Manual build

To build only the C++ components:

```bash
cd src
make
```

Use `make debug` for an unoptimized debugging build or `make clean` to remove generated objects and binaries.

## Running UniAnn

Run `uniann.sh` with a FASTA file, a coding-frame score file, and a site score file:

```bash
./bin/uniann.sh \
  --fasta genome.fa \
  --coding coding.txt \
  --scores scores.txt
```

Or, if the scripts remain beside the source tree:

```bash
./scripts/uniann.sh \
  -f genome.fa \
  -c coding.txt \
  -s scores.txt
```

### Command-line options

| Option | Description |
|---|---|
| `-f`, `--fasta FILE` | FASTA file containing a single sequence |
| `-c`, `--coding FILE` | Per-base coding-frame score file |
| `-s`, `--scores FILE` | Training/annotation file with site scores |
| `-m`, `--mult NUMBER` | Probability rescaling multiplier, normally in the range `1`–`100`; default is `e` |
| `-n`, `--noviterbi` | Suppress the Viterbi matrix diagnostic output in `out.err` |
| `-v`, `--verbose` | Enable shell tracing |
| `-h`, `--help` | Display usage information |

## Input formats

### FASTA

UniAnn accepts a single-sequence FASTA file:

```text
>contig1
ATGCGTAACTGATCGTAGCATG...
```

The FASTA header's first whitespace-delimited identifier is used as the GFF3 `seqid`.

### Coding-frame file

`coding.txt` is read by `convert_coding_frame.pl` and converted into `out.ps.txt`. The converter expects tabular, per-base values and produces five emission fields plus the base:

```text
position    noncoding_score    frame0_score    frame1_score    frame2_score    intron_score    base
```

The generated `out.ps.txt` is zero-based and has this general form:

```text
0  -1.23  -4.56  -2.34  -3.21  -1.23  A
1  -1.10  -4.01  -2.18  -3.44  -1.10  T
```

Existing `out.ps.txt` files are reused; delete this file to regenerate it.

### Site score file

The `--scores` file must contain rows for plus-strand sites. UniAnn recognizes the following site types:

- `start`
- `stop`
- `donor`
- `acceptor`

Expected whitespace/tabular fields include a sequence identifier, 1-based position, strand, site type, motif, and probability:

```text
seqid  position  strand  type       motif  probability
contig1  1         +       start      ATG    0.95
contig1  40        +       donor      GT     0.80
contig1  120       +       acceptor   AG     0.90
contig1  180       +       stop       TGA    0.98
```

UniAnn currently reads a required six-column layout and also accepts an optional seventh probability column. The final probability column is used. Only rows on the `+` strand are processed by `uniann.sh`.

## Pipeline overview

Running `uniann.sh` performs the following steps:

1. Validates the FASTA, coding-frame, and site score inputs.
2. Converts coding-frame probabilities into log-space emissions with `convert_coding_frame.pl`.
3. Scales and separates start, stop, donor, and acceptor probabilities.
4. Writes temporary score files:
   - `out.ps.txt`
   - `out.gt.txt`
   - `out.ag.txt`
   - `out.atg.txt`
   - `out.stop.txt`
5. Runs `uniann` to find the optimal state path with Viterbi dynamic programming.
6. Pipes raw GFF3 output through `gffread --tlf`.
7. Filters models with fewer than `200` CDS bases when they contain one or two exons.
8. Writes the final GFF3 annotation to:

```text
<FASTA>.uniann.gff
```

For example:

```text
genome.fa.uniann.gff
```

Viterbi diagnostics are written to `out.err` by default. Use `-n` to send them to `/dev/null`.

## Direct decoder usage

The C++ decoder can be invoked directly when all intermediate score files have been created:

```bash
./bin/uniann \
  genome.fa \
  out.ps.txt \
  out.gt.txt \
  out.ag.txt \
  out.atg.txt \
  out.stop.txt \
  > genome.gff3 2> uniann.err
```

The direct output is GFF3. Internally, the decoder:

- Reads a single FASTA sequence
- Loads per-base emissions and sparse ATG/GT/AG/stop scores
- Performs safety checks that GT, AG, ATG, and stop coordinates match the sequence
- Applies frame and minimum-length constraints
- Backtraces the optimal Viterbi path
- Emits CDS and intron features, plus internal state metadata

## Splice-site model tool

`compute_markov_scores` is the C++ splice-site score utility:

```bash
./bin/compute_markov_scores genome.fa positive_model.pwm negative_model.pwm
```

It expects zoeHMM-style model tables and writes:

- `out.gt.txt` — donor/GT scores
- `out.ag.txt` — acceptor/AG scores

The implementation loads 0HMM, 1HMM, and 2HMM tables for donors and acceptors, when available, and may combine positive and negative model scores.

## Output

A typical GFF3 result starts with:

```text
##gff-version 3
contig1	UniAnn	CDS	101	180	12.34	+	.	Parent=genome.fa.1;state=E0
contig1	UniAnn	intron	181	220	11.20	+	.	Parent=genome.fa.1;state=I0
contig1	UniAnn	CDS	221	300	13.10	+	.	Parent=genome.fa.1;state=E0
```

Feature coordinates follow the decoder's state transitions. The `state` attribute records the frame or intron state; feature scores are derived from normalized Viterbi path scores.

## Creating a source distribution

From the repository root:

```bash
./make_dist.sh
```

This creates `uniann.1.0.0.tgz` from the license, installer, executables, source code, and scripts. The version is currently fixed in `make_dist.sh` as `1.0.0`.

## Troubleshooting

### No executable output

Build the project and verify the required binaries exist:

```bash
./install.sh
ls -l bin
```

### Input files are rejected

Confirm that:

- The FASTA is non-empty and contains at least one sequence
- `-c` points to a non-empty coding-frame file
- `-s` points to a non-empty site score file
- Site positions are 1-based
- GT, AG, ATG, and stop coordinates correspond to the sequence

### Unexpected or missing predictions

Check `out.err` for the Viterbi DP/backtrace matrices. The pipeline can also emit coordinate warnings when motif positions do not match the FASTA sequence.

The default result is sensitive to probability scaling and the supplied training/site scores. Try a different `-m` multiplier when its scale is incompatible with the input model.

## Development

```bash
cd src
make debug
make clean
make
```

The C++ targets use only the C++ standard library. Perl and Bash scripts orchestrate preprocessing and GFF post-processing.

## License

UniAnn is released under the **GNU General Public License v3.0**. See [LICENSE](LICENSE) for the complete terms.


