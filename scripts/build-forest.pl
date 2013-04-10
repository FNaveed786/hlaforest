#!/usr/bin/env perl

use warnings;
use strict;
use Getopt::Long;
use HLATree;
use SamReader;
use Alignment;
#use Data::Dumper;
use Storable;

my ($refFile, $readFile, $forest_prefix, $verbose, $num_alignments);
my $block_size = 10000;


GetOptions ( 
       "ref:s"=>\$refFile, 
       "reads=s"=>\$readFile, 
       "o|out|forest=s"=>\$forest_prefix, 
       "v|verbose"=>\$verbose, 
       "n|num_alignments:i" => \$num_alignments,
       "b|block_size:i"=>\$block_size,
   );


my @hlaforest;

my $samReader = SamReader->new($readFile);
my $header = $samReader->getSamHeader;
my $alignmentSetCount = 1;
my $block_count = 0;
while (my $alignmentPtr = $samReader->getNextAlignmentSet) {
## Build a forest from this set of reads
    my $readTree = HLATree->new();
    $readTree->buildTreeFromAlignmentSet($alignmentPtr);
    push (@hlaforest, $readTree);
    $alignmentSetCount++;
#print Dumper(@hlaforest);
#    print "Processed ".$alignmentSetCount."\n" if $verbose;

    if($num_alignments) {
        last if ($alignmentSetCount >= $num_alignments);
    }
    if ($alignmentSetCount % $block_size == 0) {
        my $result = eval { store( \@hlaforest, $forest_prefix.++$block_count.".forest") };
        if( $@ )
            { warn "Serious error from Storable: $@" }
        elsif( not defined $result )
            { warn "I/O error from Storable: $!" }
        undef @hlaforest;
        my @hlaforest;
    }
}

my $result = eval { store( \@hlaforest, $forest_prefix.++$block_count.".forest") };

if( $@ )
    { warn "Serious error from Storable: $@" }
elsif( not defined $result )
    { warn "I/O error from Storable: $!" }


