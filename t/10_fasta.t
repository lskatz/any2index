#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::More tests=>2;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;
use File::chdir; # lets you modify $CWD for local chdir
use File::Which qw/which/;
use FindBin qw/$RealBin/;
use lib "$RealBin/../lib/perl5";

$ENV{PATH}="$RealBin/../scripts:$ENV{PATH}";
my $asm = "$RealBin/NC001416.fasta";
my $tempdir = tempdir(basename($0).".XXXXXX", TMP=>1, CLEANUP=>1);

subtest 'executables' => sub{
  if(!which("any2index")){
    BAIL_OUT("Could not find any2index in the path");
  }

  my @executable = qw(samtools makeblastdb bowtie2-build);
  plan tests => scalar(@executable);
  for my $exe(@executable){
    my $path = which($exe);
    is(-x $path, 1, "$exe is executable");
  }
};

subtest 'fasta' => sub{

  # Set up the working temp directory
  symlink($asm, $tempdir."/".basename($asm));
  $CWD = $tempdir;
  my $asmLink = basename($asm);

  my $command = "any2index '$asmLink' > any2index.log 2>&1";
  my $exit_code = system($command);
  is($exit_code, 0, "$command");
  if($exit_code){
    note `cat any2index.log`;
  }
};

