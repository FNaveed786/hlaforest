#!/usr/bin/env perl

use Getopt::Long;
use Bio::SeqIO;
use HLATree;
use Math::Random qw(random_normal random_binomial random_uniform_integer);

my $file;
my $numHaplotypes = 2;
my $outFasta;
my $logFile;
my $debug;
my $geneList;
my $minRefLength = 200;

GetOptions ("f|file=s"=>\$file, "o|outFasta=s"=>\$outFasta, "d|debug"=> \$debug, "g|genes:s"=>\$geneList, "m|minRefLength:i"=>\$minRefLength, "l|log=s"=>\$logFile);


my $seqIO = Bio::SeqIO->new(-file=> $file, -format => "fasta");

my %seqDb; # Stores sequences as $seqDb{$geneId}{$seqName} = sequence
while (my $seqObj = $seqIO->next_seq()) {
#    $self->addRefNode($seqObj->display_id(), $seqObj->seq());
# Parse the description field
    my ($templateLength, $imgtId, $hlaPtr) = HLATree::_extractHLAfromDesc($seqObj->display_id());
    $seqDb{$hlaPtr->[0]}{$seqObj->display_id()} = $seqObj->seq();
}

open($outFastaFh, ">", "$outFasta");
open (my $logFh, ">", "$logFile");


# Pick two random haplotypes from each gene class and generate a random number of reads
my @genes;
if ($geneList) {
    @genes = split /,/, $geneList;

}
else {
    @genes = keys %seqDb;
}

foreach my $gene (@genes) {
    my $numSeqs = keys %{$seqDb{$gene}};
    foreach (1..$numHaplotypes) {

        # Pick random haplotypes
        my $randIndex = int(rand($numSeqs));
        my %seqHash = %{$seqDb{$gene}};
        my @seqNameArray= keys %seqHash;
        my $seqName = $seqNameArray[$randIndex];
        my $seq = $seqDb{$gene}{$seqName};
        my $seqLength = length($seq);

        # Check to make sure randomly selected haplotype exceeds insert mean length
        while ($seqLength < $minRefLength) {
            print "I am stuck because $gene is only $seqLength\n";
            $randIndex = int(rand($numSeqs));
            %seqHash = %{$seqDb{$gene}};
            @seqNameArray= keys %seqHash;
            $seqName = $seqNameArray[$randIndex];
            $seq = $seqDb{$gene}{$seqName};
            $seqLength = length($seq);
        }

        print $logFh "$seqName\n";

        # Print the selected haplotypes
        print $outFastaFh ">$seqName\n";
        print $outFastaFh "$seq\n";
    }
}

close ($readOneFh);
close ($logFh);
