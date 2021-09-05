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
my $bam = "$RealBin/sorted.bam";
my $tempdir = tempdir(basename($0).".XXXXXX", TMP=>1, CLEANUP=>1);

subtest 'executables' => sub{
  if(!which("any2index")){
    BAIL_OUT("Could not find any2index in the path");
  }

  my @executable = qw(samtools);
  plan tests => scalar(@executable);
  for my $exe(@executable){
    my $path = which($exe);
    is(-x $path, 1, "$exe is executable");
  }
};

subtest 'sam' => sub{

  ## Set up the working temp directory
  
  # Create a sam file
  my $samTarget = $tempdir."/".basename($bam,".bam").".sam";
  my $exit_code = system("samtools view -h '$bam' > '$samTarget' 2> '$samTarget.log'");
  if($exit_code){
    diag "samtools view log:";
    diag `cat '$samTarget.log'`;
    BAIL_OUT("Failed to make sam from bam file");
  }
  my $localSam = basename($samTarget);
  $CWD = $tempdir;

  my $log = "any2index.log";
  my $command = "any2index '$localSam' >$log 2>$log";
  $exit_code = system($command);
  is($exit_code, 0, "$command");
  if($exit_code){
    diag "any2index log:";
    diag `cat any2index.log`;
    BAIL_OUT("Failed to run any2index");
  }
  note `ls -l .`;
  note `cat any2index.log`;
};

