package SamReader;

use Alignment;
use constant QUAL33 => 1;
use constant QUAL64 => 2;

sub new {
    my $type = shift;
    my $self = {};
    my $samFile;
    my $fh;

    if (@_) {
        $samFile = shift;
        open ($fh, "<", $samFile);
    }
    else {
        open ($fh, "<-");
    }
    $self->{FH} = $fh;
    $self->{NEXTALIGNMENT} = undef;
    $self->{SAMHEADER} = undef;


    bless($self, $type);
    return $self;
}

sub nextAlignment { 
    my $self = shift;
    my $alignmentBuffer = $self->_alignmentBuffer;
    if ($alignmentBuffer) {
        $self->_clearAlignmentBuffer;
        return $alignmentBuffer;
    }
    else {
        if (my $nextLine = $self->_nextLine) {
        my $alignment = Alignment->new($nextLine);
#        $alignment->calculateSMMQ;
        return $alignment;
        }
        else {
            return undef;
        }
    }
}

sub _clearAlignmentBuffer {
    my $self = shift;
    $self->{NEXTALIGNMENT} = undef;
}

sub _alignmentBuffer {
    my $self = shift;
    if (@_) {
        $self->{NEXTALIGNMENT} = shift;
    }
    return $self->{NEXTALIGNMENT};
#    if ($self{NEXTALIGNMENT}) {
#        return $self->{NEXTALIGNMENT};
#        $self{NEXTALIGNMENT} = $self->{SEQIO}->next_seq();
#    }
#    else {
#        $self{NEXTALIGNMENT} = $self->{SEQIO}->next_seq();
#        return $self->{NEXTALIGNMENT};
#    }
}

sub _nextLine {
    my $self = shift;
    my $fh = $self->{FH};
    my $line = <$fh>;
    return($line);
}

### This should be replaced by Alignment::qname
sub _readNameFromSam {
#    my $self = shift;
    my $alignment = shift;
#    print "DEBUG _readNameFromSam raw data: $alignment\n";
    my @split = split "\t", $alignment;
#    print "DEBUG _readNameFromSam: read name is ". $split[0] ."\n";
    return $split[0];
}
###

sub getSamHeader {
    my $self = shift;
#    print "DEBUG: getSamHeader called\n";
    while (my $line = $self->_nextLine) {
#    print "DEBUG: getSamHeader read $line \n";
        if ($line=~/^@/) {
            $self->{SAMHEADER}.=$_;
        }
        else {
            $self->_alignmentBuffer(Alignment->new($line));
#            $self->_alignmentBuffer->calculateSMMQ;
            last;
        }
    }
    return $self->{SAMHEADER};
}

sub getNextAlignmentSet {
    my $self = shift;
    my @alignments = ();
#    my $firstAlignmentInSet;

    if ( my $firstAlignmentInSet = $self->nextAlignment )  {
        push(@alignments, $firstAlignmentInSet);
        while (my $thisAlignment = $self->nextAlignment) {
            if ($thisAlignment->qname eq $firstAlignmentInSet->qname) {
                push(@alignments, $thisAlignment);
            }
            else {
                $self->_alignmentBuffer($thisAlignment);
                last;
            }
        }
    }

    if (scalar @alignments) {
    return \@alignments;
    }
    else {
        print "DEBUG in SamReader.pm getNextAlignmentSet: No more alignments found, returning undef\n";
        return undef;
    }
}

sub setQual33{ 
    my $self = shift;
    my $header = $self->getSamHeader;
    my ($convertCount, $total, $skipped) = (0, 0, 0);

    while (my $alignment = $self->nextAlignment) {
        my $qualType = $alignment->checkQualType;
        if($qualType == QUAL33) {
            $alignment->print;
        }
        elsif ($qualType == QUAL64) {
            $alignment->convertQual64ToQual33;
            $alignment->print;
            $convertCount++;
        }
        else {
            if (length($alignment->seq) > 50 ) {
                $alignment->print;
                print STDERR "Quality values ambiguous, but read length > 50 so assuming this is from the hiseq. (This only works for vaccination data)\n";
            }
            else {
                print STDERR "Excluding Alignment: " . $alignment->toString ."\n";
                $skipped++;
            }
        }
        $total++;
    }
    print STDERR "Out of $total alignments, $convertCount were converted and $skipped were discarded due to quality ambiguities\n";
}

return(1);
