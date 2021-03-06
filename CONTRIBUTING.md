# Contributing

Thank you for contributing!
This is a collaborative environment.

To contribute, fork this repo into your space and then make changes into a new branch.

## New indexer

The easiest thing to contribute is a new indexer, by adding onto the `@indexes` array at the top of the script:

    my @indexes = (
      {
        name    => "samtoolsFaidx",
        format  => "fasta",
        command => "samtools faidx __FILE__",
        help    => "Indexes for the samtools package",
        querying=> "samtools faidx __FILE__ contig:start-stop",
        version => "samtools --version-only | sed 's/^/Samtools /'",
        citation=> "echo 'Samtools citation at http://www.htslib.org'",
      },

In this example, we have a method we have called `samtoolsFaidx`
whose command is `samtools faidx __FILE__`.
`__FILE__` is evaluated to the input file specified on the command line.
This method works on files whose format is `fasta`.
A description is in the `help` key.
A tip on how to query the new index is shown in `querying`.

In your contribution, be sure to think about what the best defaults are.
Not all bioinformatics tools have baseline sane parameters.

## more indices

These are some ideas for more indexers

| format | index |
|--------|-------|
| fasta  | bwa build |
| fasta  | smalt index |
| bed | bedtools |
| vcf.gz | bcftools index |


## more ideas

* flag to ignore files it can't index
* use sane defaults for all methods

