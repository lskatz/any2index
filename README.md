# any2index

Index any file for bioinformatics.
Contributions welcome.

# Installation

Set your `PATH` to include the scripts directory

# Usage

    any2index *.fasta
    any2index --help

# Indexes available

| format | index |
|--------|-------|
| fasta  | samtools faidx |
| fasta  | bowie2-build |

# contributions welcome

## more indices

| format | index |
|--------|-------|
| fasta  | bwa build |
| fasta  | formatblastdb |
| fasta  | smalt index |
| bed | bedtools |
| sam  | samtools view -b && sort && index |
| bam  | samtools sort && index |
| vcf.gz | bcftools index |

## more ideas

* temporary folder with all indices
Strategy is to run this script in the background while using the temp dir. 
Ctrl-c to clean up the folder. 
* flag to ignore files it can't index
* flag to specify exact index to run
  * secondary flag for in depth vs fast
  * tertiary flag to bring in custom options.
    Or some kind of config. 
* use sane defaults for all methods
