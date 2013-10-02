package HLANode;
use base Node;
use Alignment;
use List::MoreUtils qw(pairwise);
use constant INFINITY => 999999999;

sub new {
    my $type = shift;
    my $self = $type->SUPER::new();
    $self->{SEQ}=undef;
    $self->{ID}=undef;
    $self->{SMMQ} = undef;
    $self->{WEIGHT} = undef;
    $self->{COVERAGE} = undef;
    $self->{READWEIGHTS} = undef;
    $self->{READWEIGHTARRAY} = undef;
    $self->{ALIGNMENTS} = undef;

    bless ($self,$type);
    return $self;
}

# Given two hla nodes, adds the weight of children to self node
sub addChildrenWeights {
    my $self = shift;
    my $add_node = shift;

    my $weight = $self->weight() + $add_node->weight();
    $self->weight($weight);

    my @children;
    my $children_ptr = $add_node->children;

    if ($children_ptr) {
        @children = @$children_ptr;

        # for each child of the specified tree's root node
        foreach my $child_ptr (@children) {
            my $this_child = $self->getChild($child_ptr->id());

            if($this_child) {
                $this_child->addChildrenWeights($child_ptr);
            }

            else {
                my $new_child = HLANode->new();
                $new_child->id($child_ptr->id());
                $new_child->weight($child_ptr->weight());
                $new_child->parent($self);
                $self->addChild($new_child);
                
                $new_child->addChildrenWeights($child_ptr);
            }
        }
    }
}

# Takes in refname (and optionally new smmq value) and returns the smmq value associated
# with the sequence at that node
sub smmq {
    my $self = shift;
    if (@_) {
        $self->{SMMQ} = shift;
    }
    return $self->{SMMQ};
}

sub isHLA {
    print "I am HLANode!\n";
}

### addAlignment adds an alignment object to an array of alignments
sub addAlignment {
    my $self = shift;
    my $alignment = shift;
    my @alignments = @{$self->{ALIGNMENTS}}; 
    push (@alignments, $alignment);
    $self->{ALIGNMENTS} = \@alignments;
}

## alignments returns an array pointer to all the alignments at the node

sub alignments {
    my $self = shift;
    return $self->{ALIGNMENTS};
}

sub delete_alignments {
    my $self = shift;
    $self->{ALIGNMENTS} = undef;
}


### getMinSMMQ returns the minimum SMMQ of all alignments stored within this Node
sub getMinSMMQ {
    my $self = shift;
    my @alignments = @{$self->alignments};
    my $minSMMQ;
    if (scalar(@alignments)) {
        my $minSMMQ = INFINITY;
        foreach (@alignments) {
            my $thisSMMQ = $_->smmq;
            $minSMMQ = $thisSMMQ if ($thisSMMQ < $minSMMQ);
        }
#print "DEBUG in HLANode.pm getMinSMMQ: Found smmq of $minSMMQ for node ".$self->id."\n";
        return $minSMMQ;
    }
    else {
#        print STDERR "ERROR in HLANode.pm: Attempting to get SMMQ values of node ".$self->id.", but it has no associated Alignments.\n";
        return INFINITY;
    }
}

### getSMMQSum returns the sum of all SMMQs of the alignments stored within this Node
sub getSMMQSum {
    my $self = shift;
    my @alignments = @{$self->alignments};
    my $minSMMQ;
    if (scalar(@alignments) > 2) {
        print STDERR "An alignment should not align to the same reference sequence twice\n";
        exit(0);
    }
    if (scalar(@alignments)) {
        my $SMMQ_sum = 0;
        foreach (@alignments) {
            my $thisSMMQ = $_->smmq;
            $SMMQ_sum += $thisSMMQ;
        }
#print "DEBUG in HLANode.pm getMinSMMQ: Found smmq of $minSMMQ for node ".$self->id."\n";
        return $SMMQ_sum;
    }
    else {
#        print STDERR "ERROR in HLANode.pm: Attempting to get SMMQ values of node ".$self->id.", but it has no associated Alignments.\n";
        return INFINITY;
    }
}



### weight returns a hash pointer to all the weightuences
sub weight {
    my $self = shift;
    if (@_) {
        $self->{WEIGHT} = shift;
##    print "DEBUG in HLANode.pm weight: Setting weight for ".$self->lineage." to ".$self->{WEIGHT}."\n";
    }
    return $self->{WEIGHT};
}

# Adds the weight value of a read to an array associated with the sequence
sub addReadWeight {
    my $self = shift;
    my $refName = shift;
    my $readWeight = shift;
    my @nodeWeights;
    if (my $nodeWeightPtr = $self->{READWEIGHTS}->{$refName}) {
        @nodeWeights = @$nodeWeightPtr;
    }
    push (@nodeWeights, $readWeight);
    $self->{READWEIGHTS}->{$refName} = \@nodeWeights;
}

sub getReadWeight {
    my $self = shift;
    my $readWeightName = shift;
    my $nodeReadWeightPtr = $self->{COVERAGE};

    unless ($readWeightName) {
        print STDERR "Error in HLANode.pm: Requested a readWeight without supplying a readWeight name\n";
        exit(0);
    }

    if (my $readWeight = $nodeReadWeightPtr->{$readWeightName}) {
        return $readWeight;
    }
    else {
        print STDERR "ERROR in HLANode.pm getReadWeight: Requested a readWeight with readWeightName $readWeightName that could not be found\n";
    }
}

sub readWeights {
    my $self = shift;
    if (@_) {
        $self->{READWEIGHTS} = shift;
    }
    return $self->{READWEIGHTS};
}

sub readWeightArray {
    my $self = shift;
    if (@_) {
        $self->{READWEIGHTARRAY} = shift;
    }
    return $self->{READWEIGHTARRAY};
}



### addWeight adds a weight to the hash of weights available to that node
# 
sub addWeight {
    my $self = shift;
    my $weightName = shift;
    my $weight = shift;
    my %nodeWeights;
    if (my $nodeWeightPtr = $self->{SEQ}) {
        %nodeWeights = %$nodeWeightPtr;
    }
    $nodeWeights{$weightName} = $weight;
    $self->{SEQ} = \%nodeWeights;
}

sub getWeight {
    my $self = shift;
    my $weightName = shift;
    my $nodeWeightPtr = $self->{SEQ};

    unless ($weightName) {
        print STDERR "Error in HLANode.pm: Requested a weight without supplying a weight name\n";
        exit(0);
    }
    if (my $weight = $nodeWeightPtr->{$weightName}) {
        return $weight;
    }
    else {
        print STDERR "ERROR in HLANode.pm: Requested a weight with weightName $weightName that could not be found\n";
    }
}

### initCoverage creates an array of 0's of length ($self->SEQ->{seqName})
sub initCoverage {
    my $self = shift;
    my $seqName = shift;
    my @coverage = ((0) x length($self->getSeq($seqName)));
    $self->{COVERAGE}->{$seqName} = \@coverage;
#    print "DEBUG in HLANode.pm initCoverage: Initializing ".$self->lineage." to\n".join(" ", @coverage)."\n";

#    if (exists $self->{COVERAGE}->{$seqName}) { print "DEBUG in HLANode.pm initCoverage: Initializing ".$self->lineage." to\n".join(" ", @{$self->{COVERAGE}->{$seqName}})."\n";}
}

### coverage returns a hash pointer to all the coverages
sub coverage {
    my $self = shift;
    return $self->{COVERAGE};
}

sub addCoverage {
    my $self = shift;
    my $coverageName = shift;
    my $coverageArrayPtr = shift;
    my @coverageArray = @{$coverageArrayPtr};
    #if (my $selfCovereagePtr = $self->{COVERAGE}->{$coverageName}) {
    if ( $self->getCoverage($coverageName)) {
        #my @selfCoverage = @{$selfCoveragePtr};
        my @selfCoverage = @{$self->{COVERAGE}->{$coverageName}};
        @selfCoverage = pairwise {$a + $b} @selfCoverage, @coverageArray;
        $self->{COVERAGE}->{$coverageName} = \@selfCoverage;
    }
    else { 
        my $debugPtr = $self->getCoverage($coverageName);
        die "Error in HLANode.pm addCoverage: addCoverage called on a node that has not been initialized\n";
    }
}

sub getCoverage {
    my $self = shift;
    my $coverageName = shift;
    my $nodeCoveragePtr = $self->{COVERAGE};

    unless ($coverageName) {
        print STDERR "Error in HLANode.pm: Requested a coverage without supplying a coverage name\n";
        exit(0);
    }

    if (my $coverage = $nodeCoveragePtr->{$coverageName}) {
        return $coverage;
    }
    else {
        print STDERR "ERROR in HLANode.pm getCoverage: Requested a coverage with coverageName $coverageName that could not be found\n";
    }
}



### seq returns a hash pointer to all the sequences
sub seq {
    my $self = shift;

    return $self->{SEQ};

}

### addSeq adds a sequence to the hash of sequences available to that node
# 
sub addSeq {
    my $self = shift;
    my $seqName = shift;
    my $seq = shift;
    my %nodeSeqs;
    if (my $nodeSeqPtr = $self->{SEQ}) {
        %nodeSeqs = %$nodeSeqPtr;
    }
    $nodeSeqs{$seqName} = $seq;
    $self->{SEQ} = \%nodeSeqs;
}

# Given the name of a sequence, return that sequence from the node
sub getSeq {
    my $self = shift;
    my $seqName = shift;
    my $nodeSeqPtr = $self->{SEQ};

    unless ($seqName) {
        print STDERR "Error in HLANode.pm: Requested a sequence without supplying a sequence name\n";
        exit(0);
    }
    if (my $seq = $nodeSeqPtr->{$seqName}) {
        return $seq;
    }
    else {
        print STDERR "ERROR in HLANode.pm: Requested a sequence with seqName $seqName that could not be found\n";
    }
}

### This is a deprecated method that assumes only a single sequence per
## node.
#sub seq {
#    my $self = shift;
#    if (@_) {
#        $self->{SEQUENCE}=shift;
#    }
#    return $self->{SEQUENCE};
#}

sub id {
    my $self = shift;
    if (@_) {
        $self->{ID}=shift;
    }
    return $self->{ID};
}

sub lineage {
    my $self = shift;
    my $lineageString = "";
    if (my $parent = $self->parent) {
        $lineageString = $parent->lineage.":".$self->id;
    }
    else {
        $lineageString = "ROOT";
    }
    return $lineageString;
}

### This simple test method should be deprecated. It is not general enough and _addNode
# should probably a function of HLATree rather than HLANode
#sub _addNode {
#    my $self = shift;
#    my $parent = shift;
#    my $hlaPtr = shift;
#    my $seq = shift;
#    my @hla = @$hlaPtr;
#    my $nodeName = shift @hla;
#
#    $self->id($nodeName);
#
#    $parent->_childExists($nodeName);
#
## Here check if any children of the parent have the same id as 
## the current id
#
## If it does, return that node
#    if (@hla) {
#        $self->_addNode($self,\@hla,$seq);
#    }
## else, return null
#    else {
#        $self->seq($seq);
#        print "adding seq $seq\n";
#    }
#}
#
#sub _childExists {
#    my $self=shift;
#    my $name=shift;
#    my @children = @$self->children();
#    foreach my $child (@children) {
#        return $child if ($child->id() eq $name);
#    }
#}

sub getChild {
    my $self=shift;
    my $name=shift;
    my @children = @{$self->children()};
    my $returnable;

#    print STDERR "HLANode.pm getChild: @children\n";

    foreach my $child (@children) {
        $returnable = $child if ($child->id() eq $name);
    }
    if ($returnable) {
        return $returnable
    }
    else {
        return undef;
    }
}


# given a single id, remove that child
sub removeChild {
    my $self = shift;
    if (@_) {
        my $child_id_to_remove = shift;
        my $childrenPtr = $self->{CHILDREN};
#        for (my $i = 0; $i < scalar(@{$childrenPtr}); $i++) {
#            if ($chidrenPtr->[$i]->id eq $child_id_to_remove) {
#                delete $childrenPtr->[$i];
#            }
#        }
        my @newChildren;
        foreach my $child (@$childrenPtr) {
            unless ($child->id eq $child_id_to_remove) {
                push @newChildren, $child;
            }
        }
        $self->{CHILDREN} = \@newChildren;
    }
    else {
        print STDERR "ERROR: removeChild called in Node.pm without an argument\n";
    }
}

# Given a single id, remove all but that child
sub removeAllButChild {
    my $self = shift;
    if (@_) {
        my $child_id_to_keep = shift;
        my @newChildren;
        push @newChildren, $self->getChild($child_id_to_keep);
        $self->{CHILDREN} = \@newChildren;
#        my $childrenPtr = $self->{CHILDREN};
#        for (my $i = 0; $i < scalar(@{$childrenPtr}); $i++) {
#            unless ($childrenPtr->[$i]->id eq $child_id_to_keep) {
#                print STDERR "Deleting ".$childrenPtr->[$i]->id." from ". $self->id."\n";
#                delete $childrenPtr->[$i];
#            }
#        }
    }
    else {
        print STDERR "ERROR: removeChild called in Node.pm without an argument\n";
    }
}
# Given an array of ids, remove all but those children
sub removeAllButChildren {
    my $self = shift;
    if (@_) {
        my $children_ids_to_keep_ptr = shift;
        my @children_ids_to_keep = @{$children_ids_to_keep_ptr};
        my @newChildren;
        foreach my $child_id_to_keep (@children_ids_to_keep) {
            push @newChildren, $self->getChild($child_id_to_keep);
        }
        $self->{CHILDREN} = \@newChildren;
    }
    else {
        print STDERR "ERROR: removeChildren called in Node.pm without an argument\n";
    }
}




return(1);

