#!/usr/bin/env perl

use warnings;
use strict;
use HLATree;
use Storable qw(retrieve);
use Getopt::Long;

my $treeFile;
my $sumFlag;
my $covFlag;
my $readWeightsFlag;
GetOptions("t=s"=>\$treeFile, "sum"=>\$sumFlag, "cov"=>\$covFlag, "readWeights"=>\$readWeightsFlag);

my $hlatree = eval { retrieve( $treeFile ) };

if ($sumFlag) {
$hlatree->printCoverageSums($hlatree->root);
}
elsif ($covFlag) {

$hlatree->printCoverages($hlatree->root);
}
elsif($readWeightsFlag) {

$hlatree->printReadWeightArray($hlatree->root);
}
