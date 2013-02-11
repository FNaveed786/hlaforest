#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use HLATree;
use SamReader;
use Alignment;
#use Data::Dumper::Simple;
use Storable;

my ($refFile, $readFile, $treeFile, $verbose);

GetOptions ("ref=s"=>\$refFile, "reads=s"=>\$readFile, "o|out|tree=s"=>\$treeFile, "v|verbose"=>\$verbose);


my $hlatree = HLATree->new();
$hlatree->buildTreeFromImgtDB_nospace($refFile);

print "ref tree built\n";

my $samReader = SamReader->new($readFile);
my $header = $samReader->getSamHeader;
my $alignmentSetCount = 1;
while (my $alignmentPtr = $samReader->getNextAlignmentSet) {
## Build a tree from this set of reads
    my $readTree = HLATree->new();
    $readTree->buildTreeFromAlignmentSet($alignmentPtr);
    $hlatree->distributeWeightedCoverage($readTree);
    print "Processed ".$alignmentSetCount++."\n" if $verbose;
}
#open (my $treeFh, ">", $treeFile);
#print $treeFh Dumper($hlatree);
#close($treeFh);

my $result = eval { store( $hlatree, $treeFile ) };

if( $@ )
        { warn "Serious error from Storable: $@" }
        elsif( not defined $result )
            { warn "I/O error from Storable: $!" }


