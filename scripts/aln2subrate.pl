#!/usr/bin/perl
# Given an aln file from stdin, returns the substitution rate
# Usage:
#   aln2subrate.pl < file.aln
#
# Outputs a tab delimited file where the first value is the average substituion rate over the entire read and the remaining is the subsitution rate at that position


use strict;
use warnings;

# read through header

my $buffer = <>;

sub read_header {
    my $header;
    if ($buffer =~ /^[#@]/) {
        $header = $buffer;
        while (<>) {
            $buffer = $_;
            if ($buffer =~/^[#@]/) {
                $header .= $buffer;
            }
            else { 
                last; 
            }

        }
    }
    return $header;
}

sub get_next_aln {
    my ($id, $ref, $sim);
    if ($buffer) {
#        if ($buffer =~ /^>/) {
            $id = $buffer;
            $ref = <>;
            $sim = <>;
            $buffer = <>;
#        }
        return ($id, $ref, $sim);
    }
    else {
        return (0);
    }
}

# Removes leading and trailing whitespace
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}

sub mism_pos {
    my ($str1, $str2) = @_;
    my @mism_pos;

    for my $i (0 .. length($str1) - 1) {
        if (substr($str1, $i, 1) ne substr($str2, $i, 1) ) {
            push @mism_pos, $i;
        }
    }
    return \@mism_pos;
}

my $num_reads = 0;
my $num_errors = 0;
my $total_bases = 0;
my @errors;
my @error_rate;

my $header = read_header();

# populate error array
while (my ($id, $ref, $sim) = get_next_aln()) {
    last unless $id;
    $num_reads++;

    $ref = trim($ref);
    $sim = trim($sim);

    my @diffs = @{mism_pos($ref, $sim)};
    foreach my $diff (@diffs) {
        $errors[$diff]++;
        if ($diff > 200) {
            print "$ref\n";
            print "$sim\n";
        }
    }

    $total_bases += length($sim);
    $num_errors += scalar(@diffs);

}

# Calculate the error rate over each position
for( my $i = 0; $i <= $#errors; $i++) {
    if ($errors[$i]) {
        $error_rate[$i] = $errors[$i] / $num_reads;
    }
    else {
        $error_rate[$i] = 0;
    }

}

#print "inferred read length: ";
#print scalar(@errors);
#print "\n";


print $num_errors / $total_bases;
print "\t";
print join "\t", @error_rate;
print "\n";

