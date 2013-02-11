#!/usr/bin/env perl

# Given a file with interspersed fastq 33 and fastq 64, converts all  qual strings to fastq 33
use SamReader;
use Getopt::Long;
use warnings;
use strict;

my $filename;

GetOptions ("f=s"=>\$filename);

my $samReader;
# Initialize samReader
if ($filename) {
    $samReader = SamReader->new($filename);
}
else {
    $samReader = SamReader->new();
}

$samReader->setQual33;
