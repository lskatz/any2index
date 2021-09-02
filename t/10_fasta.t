#!/usr/bin/env perl
use strict;
use warnings;
use Data::Dumper;
use Test::More tests=>1;
use File::Temp qw/tempdir/;
use File::Basename qw/basename/;
use File::chdir; # lets you modify $CWD for local chdir
use FindBin qw/$RealBin/;
use lib "$FindBin/../lib/perl5";

$ENV{PATH}="$RealBin/../scripts:$ENV{PATH}";
my $asm = "$RealBin/NC001416.fasta";
my $tempdir = tempdir(basename($0).".XXXXXX", TMP=>1, CLEANUP=>1);

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

