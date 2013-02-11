#!/usr/bin/env perl
#
# aggregate_sim_scores.pl
#
# Given a directory of simulations generated via a shell script
# implemented at http://pourmand-wiki.soe.ucsc.edu/doku.php?id=hyjkim:binf:hla:perl#systematic_simulations_v2
# outputs overall correct call rate and correct call rate for each
# gene. Also generates a ranked list of incorrectly called haplotypes

use warnings;
use strict;
use Getopt::Long;

my ($in_dir);

GetOptions("i|in-dir=s"=>\$in_dir);

my %correct;
my $correct_tier=0;
my %total;
my $total_tier=0;

my %correctly_called;
my %incorrectly_called;

while (my $file = <$in_dir/*/sim-score.txt>) {
#    print "Processing $file \n";
    open (my $fh, "<", $file);
    while (<$fh>) {
        my @split = split /\t/, $_;
        my $ref = shift @split;
        my $gene = shift @split;
        my $correct = pop @split;
        my $ref_num = pop @split;

        foreach (my $tier = 1; $tier <= $ref_num; $tier++) {
            # Check if correct
            if ($correct >= $tier) {
                $correct{$gene}{$tier}++;
                $correctly_called{$ref}++;
            }
            # If it's false
            else {
                $incorrectly_called{$ref}++;

            }
                $total{$gene}{$tier}++;
        }
    }
    close ($fh);
}

foreach my $gene (sort keys %total) {
    print $gene;
    foreach my $tier(sort {$a <=> $b} keys %{$total{$gene}}) {
        my $correct_rate = $correct{$gene}{$tier} / $total{$gene}{$tier};
        print "\t";
        print $correct_rate;

    }
    print "\n";
}

foreach my $haplotype (sort {$incorrectly_called{$b} <=> $incorrectly_called{$a}}keys %incorrectly_called) {
#    print join "\t", $haplotype, $incorrectly_called{$haplotype}, "\n";
}

foreach my $haplotype (sort {$correctly_called{$b} <=> $correctly_called{$a}}keys %correctly_called) {
#    print join "\t", $haplotype, $correctly_called{$haplotype}, "\n";
}
