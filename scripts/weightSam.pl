#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use HLATree;
use SamReader;
use Alignment;
#use Data::Dumper;

my ($readFile, $verbose);
GetOptions ("reads=s"=>\$readFile, "v|verbose"=>\$verbose);

my $samReader = SamReader->new($readFile);
my $header = $samReader->getSamHeader;
my $alignmentSetCount = 1;

while (my $alignmentPtr = $samReader->getNextAlignmentSet) {
    my $readTree = HLATree->new();
    $readTree->buildTreeFromAlignmentSet($alignmentPtr);
    #$readTree->printNodeWeights($readTree->root);
    my @alignments = @{$readTree->_returnAlignments($readTree->root)};
    foreach my $alignment (@alignments) {
#        print join "\t", $alignment->rname, $alignment->weight, "\n";
        $alignment->setOptionalFieldValue('XW', 'f',  $alignment->weight);
        print $alignment->toString();
    }
    print STDERR "Processed ".$alignmentSetCount++."\n" if $verbose;
}
# Take in a sam file
# Get set of alignments
# Calculate read weights
# Collate read weights with alignments in alignment set
# Print out alignment set
