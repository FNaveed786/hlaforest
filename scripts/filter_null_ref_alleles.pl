#!/usr/bin/env perl
# This simple script only prints out HLA haplotypes that do not have a null
# allele (or as defined by the hla nomenclature to contain an 'N' in the 4th
# digit)

use warnings;
use strict;
use Bio::SeqIO;
use Getopt::Long;

my $filename;
GetOptions("f|filename=s"=>\$filename);

my $seqio_object = Bio::SeqIO->new(-file => $filename);
my $seqio_out = Bio::SeqIO->new(-format => "fasta");

while (my $seq_object   = $seqio_object->next_seq) {
    my (@ref_split) = split "_", $seq_object->id;
    $seqio_out->write_seq($seq_object) unless $ref_split[1] =~ /N$/;
}
