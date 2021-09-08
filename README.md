# any2index

Index any file for bioinformatics.
Sane defaults.

# Installation

Set your `PATH` to include the scripts directory

# Usage

    any2index *.fasta *.bam
    any2index --help

To list available indexes:

    any2index --list

To get help on any one index:

    any2index --help-with bowtie2Build

# Highlighted indexes

This list might not be comprehensive

| format | index |
|--------|-------|
| fasta  | samtools faidx |
| fasta  | formatblastdb |
| fasta  | minimap2 ONT |
| fasta  | minimap2 PacBio |
| fasta  | bowie2-build |
| fasta  | bwa |
| fasta  | formatblastdb |
| fasta  | mash |
| fastq  | mash |
| bam  | samtools index |
| sam  | samtools view -b && sort && index |

# contributions welcome

_More information_: [CONTRIBUTING.md](CONTRIBUTING.md)

