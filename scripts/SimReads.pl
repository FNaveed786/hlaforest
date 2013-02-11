#!/usr/bin/env perl

use Getopt::Long;
use Bio::SeqIO;
use HLATree;
use Math::Random qw(random_normal random_binomial random_uniform_integer);

my $file;
my $numReads = 1000;
my $readLength = 50;
my $paired;
my $insertMean= 200;
my $insertSD = 50;
my $numHaplotypes = 2;
my $errorRate = 0;
my $outPrefix;
my $debug;
my $geneList;

my %errorTable = (
    "A" => ["C","G","T"],
    "C" => ["A","G","T"],
    "G" => ["C","A","T"],
    "T" => ["C","G","A"],
);

GetOptions ("f|file=s"=>\$file, "n|numReads:i"=>\$numReads, "l|length:i"=>\$readLength, "pe|paired"=>\$paired, "o|outPrefix=s"=>\$outPrefix, "i|insertMean:i"=>\$insertMean, "s|insertSD:i"=>\$insertSD, "d|debug"=> \$debug, "g|genes:s"=>\$geneList, "e|errorRate:f"=>\$errorRate);


my $seqIO = Bio::SeqIO->new(-file=> $file, -format => "fasta");

my %seqDb; # Stores sequences as $seqDb{$geneId}{$seqName} = sequence
while (my $seqObj = $seqIO->next_seq()) {
#    $self->addRefNode($seqObj->display_id(), $seqObj->seq());
# Parse the description field
    my ($templateLength, $imgtId, $hlaPtr) = HLATree::_extractHLAfromDesc($seqObj->display_id());
    $seqDb{$hlaPtr->[0]}{$seqObj->display_id()} = $seqObj->seq();
}

my ($readOneFh, $readTwoFh);
if ($paired) {
    open($readOneFh, ">", "$outPrefix\_1.fa");
    open($readTwoFh, ">", "$outPrefix\_2.fa");
}
else {
    open($readOneFh, ">", "$outPrefix.fa");
}

open (my $logFh, ">", "$outPrefix\_chosen_haplotypes.txt");


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

        my $randIndex = int(rand($numSeqs));
        my %seqHash = %{$seqDb{$gene}};
        my @seqNameArray= keys %seqHash;
        my $seqName = $seqNameArray[$randIndex];
        my $seq = $seqDb{$gene}{$seqName};
        my $seqLength = length($seq);

        while ($seqLength < $insertMean) {
            print "I am stuck because $gene is only $seqLength\n";
            $randIndex = int(rand($numSeqs));
            %seqHash = %{$seqDb{$gene}};
            @seqNameArray= keys %seqHash;
            $seqName = $seqNameArray[$randIndex];
            $seq = $seqDb{$gene}{$seqName};
            $seqLength = length($seq);
        }

        print $logFh "$seqName\n";

        for (my $i = 1; $i <=  $numReads; $i++) {
            if ($paired) {
                my $insertSize = random_normal(1, $insertMean, $insertSD);
                my $randSeqStart = int(rand($seqLength));
                if ($randSeqStart > $seqLength - $insertSize) {
                    $i--;
                    next;
                }
                my $simFrag= substr $seq, $randSeqStart, $insertSize;
                # This was before adding the error function
#                my $readOne = substr $simFrag, 0,$readLength;
#                my $readTwo = revdnacomp(substr $simFrag, -$readLength);
                my $readOne = insertErrors(substr ($simFrag, 0,$readLength), $errorRate, $readLength);
                my $readTwo = insertErrors(revdnacomp(substr ($simFrag, -$readLength)), $errorRate, $readLength);

                if ($debug) {
                    print "DEBUG R: $simFrag\n" ;
                    print "DEBUG 1: $readOne\n";
                    print "DEBUG 2: $readTwo\n\n";
                }

                print $readOneFh ">$seqName|$i|$randSeqStart\n$readOne\n";
#                print $readTwoFh ">$seqName|$i|".($randSeqStart+$insertSize-$readLength)."\n$readTwo\n";
                print $readTwoFh ">$seqName|$i|$randSeqStart\n$readTwo\n";
            }

            else {
                my $randSeqStart = int(rand($seqLength));
# Check to make sure the simulated read will correctly fall within the boundaries of the read
                if ($randSeqStart > $seqLength - $readLength) {
                    $i--;
                    next;
                }
                my $simRead = substr $seq, $randSeqStart, $readLength;
                print $readOneFh ">$seqName|$i|$randSeqStart\n$simRead\n";
            }
        }
    }
}

sub revdnacomp {
# my $dna = @_;  
# the above means $dna gets the number of 
# arguments in @_, since it's a scalar context!

    my $dna = shift; # or   my $dna = shift @_;
# ah, scalar context of scalar gives expected results.
# my ($dna) = @_; # would work, too

    my $revcomp = reverse($dna);
    $revcomp =~ tr/ACGTacgt/TGCAtgca/;
    return $revcomp;
}

## insertErrors
# Given a sequence read and error rate, randomly distribute errors throughout the read
# and return the sequence

sub insertErrors {
    my $sequence = shift;
    my $errorRate = shift;
    my $readLength = shift;

    # If the errorRate is zero, just return the sequence
    if ($errorRate > 0) {
        my $numErrors;
        my @errorPositions;

        # If readLength is not provided by caller, generate it from the sequence
        $readLength = length ($sequence) unless ($readLength);

        # generate the number of errors to be distributed
        $numErrors = random_binomial(1, $readLength, $errorRate);

        # get positions to replace
        @errorPositions = random_uniform_integer($numErrors, 0, $readLength-1);
        foreach my $errorIndex (@errorPositions) {
            my $base = substr($sequence, $errorIndex, 1);
            substr($sequence, $errorIndex, 1) = misreadBase($base);
        }
    }
    return $sequence;
}


# Given a base, return a random base that is not the same base that was given ## yay circular confusing documentation
sub misreadBase{
    my $base = shift;
    return $errorTable{$base}->[random_uniform_integer(1,0,2)];
}
