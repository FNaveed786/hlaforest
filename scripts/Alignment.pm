package Alignment;

# Set some constant values
use constant QUAL33 => 1;
use constant QUAL64 => 2;
use constant QUAL33OFFSET => 33;
use constant INDELQUAL => 20;

# This object stores the relevant objects of a sam line alignment

sub new {
    my $type = shift;
    my $samLine = shift;
    my $self = {};
    my %optional;
    

    if ($samLine) {
    chomp ($samLine);
    my @samSplit = split /\t/, $samLine;

    $self->{QNAME} = shift @samSplit;
    $self->{FLAG} = shift @samSplit;
    $self->{RNAME} = shift @samSplit;
    $self->{POS} = shift @samSplit;
    $self->{MAPQ} = shift @samSplit;
    $self->{CIGAR} = shift @samSplit;
    $self->{RNEXT} = shift @samSplit;
    $self->{PNEXT} = shift @samSplit;
    $self->{TLEN} = shift @samSplit;
    $self->{SEQ} = shift @samSplit;
    $self->{QUAL} = shift @samSplit;

## This is a minimal set of alignments needed for forest building
#    $self->{QNAME} = shift @samSplit;
#    shift @samSplit;
#    $self->{RNAME} = shift @samSplit;
#    shift @samSplit;
#    shift @samSplit;
#    $self->{CIGAR} = shift @samSplit;
#    shift @samSplit;
#    shift @samSplit;
#    shift @samSplit;
#    $self->{SEQ} = shift @samSplit;
#    $self->{QUAL} = shift @samSplit;

    foreach my $field (@samSplit) {
        my ($tag, $type, $value) = split /:/, $field;
        $optional{$tag}{TYPE} = $type;
        $optional{$tag}{VALUE} = $value;
    }
    $self->{OPTIONAL} = \%optional;

    $self->{SMMQ} = undef;
    #$self->{WEIGHT} = undef;

    bless($self, $type);
    return ($self);
    }
    else {return 0;}
}


sub seq {
    my $self = shift;
    return $self->{SEQ};
}

sub qual {
    my $self = shift;
    if (@_) {
        $self->{QUAL} = shift;
    }
    return $self->{QUAL};
}


sub qname {
    my $self = shift;
    return $self->{QNAME};
}

sub rname {
    my $self = shift;
    return $self->{RNAME};
}

sub pos {
    my $self = shift;
    return $self->{POS};
}

sub cigar {
    my $self = shift;
    return $self->{CIGAR};
}

sub getOptionalFieldValue {
    my $self = shift;
    my $optionalTag = shift;
    return $self->{OPTIONAL}->{$optionalTag}->{VALUE};
}

sub setOptionalFieldValue {
    my $self = shift;
    my $optionalTag = shift;
    my $optionalType = shift;
    my $optionalValue = shift;
    $self->{OPTIONAL}->{$optionalTag}->{TYPE} = $optionalType;
    $self->{OPTIONAL}->{$optionalTag}->{VALUE} = $optionalValue;
}

# Returns an array of arrays
# where each sub array has array[0] = numbases
# and array[1] = operationtype
sub cigarArray {
    my $self   = shift;
    my $cigar  = $self->cigar;
    my @result;

    while ($cigar =~ /(\d+)(\w)/g) {
        #my @cigarPortion = [$1, $2];
        my @cigarPortion = ($1, $2);
        push @result, \@cigarPortion;
    }
    return \@result;
}

# returns the indices of mismatched positions based on optional MD field
sub getMismatchPositions {
    my $self = shift;
    my $currentPos = 0;
    my @mismatchPositions;

    my $mdString = $self->getOptionalFieldValue("MD");
## For some reason this block wasn't working
#    while ($mdString =~ /^(\d+|[ATGC]+|\^[ATGC]+)/g) {
#        print "$1\n";
#
#    }
##
    $mdString =~ s/(\d+|[ATGC]+|\^[ATGC]+)/$1,/g;
    my @mdArray = split /,/, $mdString;
    foreach my $subMD(@mdArray) {
# If mismatch detected
        if ($subMD =~ /^[ATGC]/) {
            push @mismatchPositions, $currentPos++;
        }
# Do not add to the mismatch position here because indels are accounted for
# by the cigar string. Just modify the currentPosition by adding the number of ATGC's

# If Deletion is detected
        elsif ($subMD =~ /^\^/) {
            $currentPos += ($subMD =~ tr/ATGC//);
        }
# Correctly matched base case
        elsif ($subMD =~ /^\d/) {
            $currentPos += $subMD;
        }
# Just a crappy error catching message
        else {
            print STDERR "WARNING in Alignment.pm: Unusual type found in MD field parsing: $subMD\n";
            exit(0);
        }
    }
    return \@mismatchPositions;
}

sub weight{
    my $self = shift;
    if (@_) {
        $self->{WEIGHT} = shift;
##        print "DEBUG in Alignment.pm weight: Setting weight for ".$self->toString." to ".$self->{WEIGHT}."\n";
    }
    return $self->{WEIGHT};
}

# Checks if smmq has been calculate for this alignment. If it hasnt' calculates it and then returns the value
sub smmq {
    my $self = shift;
    unless  ($self->{SMMQ}) {
        $self->calculateSMMQ;
    }
    return $self->{SMMQ};
}

# Given an alignment, calculates the Sum of MisMatch Qualities (SMMQ) with
# indels contributing CONSTANT INDELQUAL. 
# Currently only checks for cigar types M,I and D
sub calculateSMMQ {
    my $self = shift;
    my $indelCount = 0;
    my $smmq=0;

# Count Indels from cigar string
    my @cigarArray = @{$self->cigarArray};
    foreach $subArrayPtr (@cigarArray) {
        if ($subArrayPtr->[1] eq "I" || $subArrayPtr->[1] eq "D") {
            $indelCount+=$subArrayPtr->[0];
        }
        elsif (! $subArrayPtr->[1] eq "M") {
            print STDERR "Warning in Alignment.pm calculateSMMQ: Found an unexpected cigar type ". $subArrayPtr->[1].". Ignoring\n";

        }
    }
    $smmq += INDELQUAL * $indelCount;

# Sum qualities of mismatching bases using MD:Z field
    my @mismatchPositions = @{$self->getMismatchPositions};
    my @quals = @{$self->qualArray};
    foreach my $pos (@mismatchPositions) {
        $smmq += $quals[$pos]; 
    }

    #return $smmq;
    $self->{SMMQ} = $smmq;
}

## qualArray
# returns an array of quality scores in int form
sub qualArray{
    my $self = shift;
    my @quals = split //, $self->qual;
    my @intQuals;
    foreach my $qual (@quals) {
        push @intQuals, ord($qual) - QUAL33OFFSET; 
    }
    return \@intQuals;
}

## qualArrayToString
# Given an array of numerical quality scores, return a phred33 qual string


## checkQualType
# returns 1 if qual 33, 2 if qual 64 or 0 if unknown
sub checkQualType {
    my $self = shift;
# Split the quality values into an array
    my @quals = split //, $self->qual;
    my ($is33, $is64) = (0, 0) ;

# Convert each quality value to an ord.
    foreach my $qualChar (@quals) {
        my $qual = ord ($qualChar);

        if ($qual >= 33 && $qual <= 63) {
            $is33++;
        }
        if ($qual >= 75) {
            $is64++;
        }
    }
    if ($is33 && $is64) {
        return 0;
    }
    elsif($is33) {
        return QUAL33;
    }
    elsif($is64) {
        return QUAL64;
    }
    else {
        print STDERR "Warning in Alignment.pm checkQualType: Could not distinguish between qual 33 and qual 64 for line: ";
        $self->print;
        print STDERR "Found $is33 qual 33 values and $is64 qual 64 values\n";
        return 0;
    }
}

sub convertQual33ToQual64 {
    my $self = shift;

    my @quals = split //, $self->qual;
       
    for (my $i = 0; $i <= $#quals; $i++) {
        my $qual = ord ($quals[$i]);
        $quals[$i] = chr($qual+31);
    }
    $self->qual(join "", @quals);
}

sub convertQual64ToQual33 {
    my $self = shift;
    my @quals = split //, $self->qual;
       
    for (my $i = 0; $i <= $#quals; $i++) {
        my $qual = ord ($quals[$i]);
        $quals[$i] = chr($qual-31);
    }
    $self->qual(join "", @quals);
}

### getCoverageArray
#  Returns a weighted coverage array relative to the reference sequence.
# This means deletions are given 0 coverage at that position in the reference
# and insertions are completely skipped
#
# **NOTE** The array returned by this function does not start at index 0.
# This is handled correctly by List::Util "pairwise", but if a custom
# method is used, this should be taken into account
sub getCoverageArray {
    my $self = shift;
    my @cigarArray = @{$self->cigarArray};
    my @coverage;
    my $readPosition= 0;
    my $refPosition = $self->pos - 1; # converting 1 based to 0 based position
# Generate a hash of mismatch positions for quick checking
    my %mismatchPositionHash = map {$_ =>1} @{$self->getMismatchPositions};
    my @quals = @{$self->qualArray};

    my $weight = $self->weight;
    unless ($weight) {
        die "ERROR in Alignment.pm getCoverageArray: No weight associated with this alignment.";
    }

# Iterate through cigar groups
    foreach $subArrayPtr (@cigarArray) {
        my @subArray = @{$subArrayPtr};
        if ($subArrayPtr->[1] eq "M" ) {
            for(my $start = $readPosition; $readPosition <= ($start + $subArrayPtr->[0]); $readPosition++) {
# Check if this match was a mismatch
                if (exists $mismatchPositionHash{$readPosition}) {
# Set the coverage to the probability that the read was an error
                    $coverage[$refPosition] = (10**(- $quals[$readPosition] / 10)) * $self->weight;
                    $refPosition++;
                }
# If it was not a mismatch
                else {
# Set the coverage to the probability that the read was not an error
                    $coverage[$refPosition] = (1 - 10**(- $quals[$readPosition] / 10)) * $self->weight;
                    $refPosition++;
                }
            }
        }
# Check if this base was an insertion in the read (Bases not present in reference)
        elsif($subArrayPtr->[1] eq "I") {
            my $numInsertions = ($readPosition + $subArrayPtr->[0]);
# This iterates the read position without changing the ref positions
# insertions in the read do not change the total coverage outputted here
            for(; $readPosition <= $numInsertions ; $readPosition++) { }
        }
# Check if this base was a deletion in the read (Bases are present in reference, but 
# not in the read)
        elsif($subArrayPtr->[1] eq "D") {
            my $numDeletions = ($readPosition + $subArrayPtr->[0]);
# This iterates the ref position without changing the read positions
            for (; $refPosition <= $numDeletions; $refPosition++) {
# Sets the coverage of this position to the probability of an indel based on the
# CONSTANT INDELQUAL
                    $coverage[$refPosition] = 1 - 10^(- INDELQUAL / 10);
            }
        }
# Softclipped bases are treated as insertions. They do not add to the
# coverage provided by the read.
        elsif($subArrayPtr->[1] eq "S") {
# This iterates the read position without changing the ref positions
# softclipped bases in the read do not change the total coverage outputted here
            my $numSubstitutions = $readPosition + $subArrayPtr->[0];
            for(; $readPosition <= $numSubstitutions ; $readPosition++) { }
        }
#
        else {
            die "ERROR in Alignment.pm getCoverageArray: Found an unexpected cigar type ".$subArrayPtr->[1] .". Dying rather ungracefully. \n";
        }
    }
    return \@coverage;
}

sub print {
    my $self = shift;
# Print required fields
    print join ("\t", 
            $self->{QNAME}, 
            $self->{FLAG}, 
            $self->{RNAME}, 
            $self->{POS}, 
            $self->{MAPQ}, 
            $self->{CIGAR}, 
            $self->{RNEXT}, 
            $self->{PNEXT}, 
            $self->{TLEN}, 
            $self->{SEQ}, 
            $self->{QUAL});

# Print optional fields
    foreach my $tag (keys %{$self->{OPTIONAL}}) {
        print "\t";
        print join(":", $tag, $self->{OPTIONAL}->{$tag}->{TYPE},
                $self->{OPTIONAL}->{$tag}->{VALUE});
    }
    print "\n";
}

sub toString{
    my $self = shift;
    my $string = "";
# Print required fields
    $string = join ("\t", 
            $self->{QNAME}, 
            $self->{FLAG}, 
            $self->{RNAME}, 
            $self->{POS}, 
            $self->{MAPQ}, 
            $self->{CIGAR}, 
            $self->{RNEXT}, 
            $self->{PNEXT}, 
            $self->{TLEN}, 
            $self->{SEQ}, 
            $self->{QUAL});

# Print optional fields
    foreach my $tag (keys %{$self->{OPTIONAL}}) {
        $string .= "\t";
        $string .= join(":", $tag, $self->{OPTIONAL}->{$tag}->{TYPE},
                $self->{OPTIONAL}->{$tag}->{VALUE});
    }
    $string .= "\n";
}


return(1);
