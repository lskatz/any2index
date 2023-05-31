#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::More tests=>1;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;
use File::chdir; # lets you modify $CWD for local chdir
use File::Which qw/which/;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";

$ENV{PATH}="$RealBin/../scripts:$ENV{PATH}";
my $asm = "$RealBin/NC001416.fasta";
my $tempdir = tempdir(basename($0).".XXXXXX", TMP=>1, CLEANUP=>1);

system("any2index $asm 2>/dev/null");

subtest 'indexes on fasta' => sub{
  plan tests => 2;

  
  my $type = `check4index $asm 2>/dev/null`;

  is(1, $type =~ /blastn/, "blastn index");
  is(1, $type =~ /faidx/,  "faidx index");
};

