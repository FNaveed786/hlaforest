package HLATree;
use HLANode;
use Bio::SeqIO;
use SamReader;
use Alignment;
use List::Util qw(min sum);

sub new {
    my $type = shift;
    my $self = {};
    my $root = HLANode->new();


    $self->{ROOT} = $root;
    $root->{ID}="ROOT";
    

    bless ($self, $type);
#    weaken($self->{ROOT});
    return $self;
}

sub root {
    my $self=shift;
    return $self->{ROOT};
}

sub buildTreeFromAlignmentSet {
    my $self = shift;
    my $alignmentsPtr = shift;

    foreach my $alignment (@$alignmentsPtr) {
        if ($alignment) {
            $self->addAlignNode($alignment);
        }
    }
    $self->_distributeSMMQ($self->root);
    $self->_calculateWeightFromSMMQ($self->root, 1);
    $self->clearAllAlignNodes($self->root);
}

sub addAlignNode {
    my $self = shift;
    my $alignment = shift;
    my ($templateLength, $imgtId, $hlaPtr) = _extractHLAfromDesc($alignment->rname);
    
    $self->_insertAlignNode($alignment, $self->root, $hlaPtr);
}

# Removes alignments from children of specified node to minimize file size
sub clearAllAlignNodes {
    my $self = shift;
    my $currentNode = shift;

    my @children = @{$currentNode->children()};

    $currentNode->delete_alignments();

    foreach my $child (@children) {
        $self->clearAllAlignNodes($child);
    }
}




sub _insertAlignNode {
    my $self = shift;
    my $alignment = shift;
    my $currentNode = shift;
    my $hlaPtr = shift;


# Check if we should traverse further down the tree before creating the sequence node
    if(scalar(@$hlaPtr)) {
        my $childId = shift @$hlaPtr;
# check if the children node exist, and if they add the next node
        if (my $child = $currentNode->getChild($childId)) {
            $self->_insertAlignNode($alignment, $child, $hlaPtr);
        }
# if the child node doesn't exist, create it
        else {
            my $child = HLANode->new();
            $child->parent($currentNode);
            $child->id($childId);
            $currentNode->addChild($child);
            $self->_insertAlignNode($alignment,$child, $hlaPtr);
        }
    }
#  If this is the deepest level we shall traverse, then
# add the coverage and SMMQ for this value
    else {
        $currentNode->addAlignment($alignment);
    }
}

## calculateReadWeights
# Given a read tree, calculates the read weights from the SMMQ's at each node
# and assigns it to the hlanode->weight of that node
sub calculateReadWeights {
    my $self = shift;
# set the SMMQ of all the nodes given the SMMQs of the alignments
#   such that the SMMQ returned from any node is the min(self->smmq, children->smmq)
    $self->_distributeSMMQ($self->root);
# Now calculate the weights
}

## _distributeSMMQ
# Used by calculateReadWeights to assign SMMQ values to nodes without alignments
sub _distributeSMMQ {
    my $self = shift;
    my $currentNode = shift;
    my @childrenSMMQ;

    if (scalar @{$currentNode->children}) {
        my @smmqs;
        push @smmqs, $currentNode->getMinSMMQ;

        foreach my $child (@{$currentNode->children}) {
            push @smmqs, $self->_distributeSMMQ($child);
        }

        $currentNode->smmq(min(@smmqs));
        return $currentNode->smmq;
    }
# If there are no children, set node->smmq to min of all alignment smmqs
    else {
        #$currentNode->smmq($currentNode->getMinSMMQ);
        $currentNode->smmq($currentNode->getSMMQSum);
        return $currentNode->smmq;
    }
}

## _calculateWeightFromSMMQ
# After SMMQ values have been distributed to nodes, calculate the read weights of each alignment
sub _calculateWeightFromSMMQ{
    my $self = shift;
    my $currentNode = shift;
    my $totalWeight = shift;

# Get the weights of al the children
    my @alignments= @{$currentNode->alignments};
    my @children = @{$currentNode->children};

# Convert SMMQ into weights
    my @localProb;
    my @objectPtr; #objectPtr holds pointers to both alignments and nodes
            

    foreach my $alignment(@alignments) {
        #$localProb{weight}{$alignment} = 10**(-$alignment->smmq/10);
        push (@localProb, 10**(-$alignment->smmq/10));
        push (@objectPtr, $alignment);
    }
    foreach my $child(@children) {
        #$localProb{$child} = 10**(-$child->smmq/10);
        push (@localProb, 10**(-$child->smmq/10));
        push (@objectPtr, $child);
    }

# Calculate the sum of all local probabilities
    my $probSum;
    foreach my $prob (@localProb) {
        $probSum += $prob;
    }


# Set the weights
    for (my $i = 0; $i < scalar(@objectPtr); $i++) {
        $objectPtr[$i]->weight($localProb[$i]/$probSum*$totalWeight);
    }

# Call this function recursively on children nodes
    foreach my $child (@children) {
        $self->_calculateWeightFromSMMQ($child, $child->weight);
    }
}

## buildTreeFromImgtDB_nospace
# Given a fasta file where all the spaces in the description fields are
# replaced by '_', build a reference tree that stores the sequences
# of each haplotype.
sub buildTreeFromImgtDB_nospace {
    my $self = shift;
    my $file = shift;
    my $seqIO = Bio::SeqIO->new(-file=> $file, -format => "fasta");

    while (my $seqObj = $seqIO->next_seq()) {
#        print "DEBUG in buildTreeFromImgtDB_nospace in HLATree.pm: Adding Ref Node".$seqObj->display_id()."\n";
        $self->addRefNode($seqObj->display_id(), $seqObj->seq());
    }
}

## addRefNode
# Given the description field of an ImgtDB_nospace file and the 
# nucleotide sequence, add a node to the tree
sub addRefNode {
    my $self = shift;
    my $seqName = shift;
    my $seq = shift;
    my ($templateLength, $imgtId, $hlaPtr) = _extractHLAfromDesc($seqName);

    $self->_insertRefNode($seqName, $seq, $self->root, $hlaPtr);

}

## _insertRefNode
# A private recursive function used to actually insert the Node
# into the tree. 
sub _insertRefNode {
    my $self = shift;
    my $seqName = shift;
    my $seq = shift;
    my $currentNode = shift;
    my $hlaPtr = shift;
    my @hla = @{$hlaPtr};


# First check if there are any more ID's in the HLAId($hlaPtr) array
    if (scalar(@$hlaPtr)) {
        my $childId = shift @$hlaPtr;
# check if any child exists with the same ID
# if it does, call insertNode with a shifted HLAID array on the child
        if (my $child = $currentNode->getChild($childId)) {
            $self->_insertRefNode($seqName, $seq, $child, $hlaPtr);
        }
# If not, first create the node, then call _insertRefNode on the 
# newly created child
        else {
            my $child = HLANode->new();
            $child->parent($currentNode);
            $child->id($childId);
            $currentNode->addChild($child);
            $self->_insertRefNode($seqName, $seq, $child, $hlaPtr);
        }
    }
# If there are no more ID's in the HLAId array, add the sequence to the current Node
    else {
        $currentNode->addSeq($seqName, $seq);
        $currentNode->initCoverage($seqName);
    }
}



## _extractHLAfromDesc
# Given a description from the hla_nuc_nospace.fasta file, return an array containing:
# (imgtId, readlength, $hlaPtr)
sub _extractHLAfromDesc {
    my $seqName = shift;
    my @split = split /[_:*]/, $seqName;

# Remove preceding '>HLA'
    shift @split;
# Remove final 'bp'
    pop @split;
    my $readLength = pop @split;
    my $imgtID = shift @split;
    my @hla = @split;
    return ($readLength, $imgtId, \@hla);
}

## distributeWeightedCoverage
# Called from a reference tree with a read tree sent in as an option.
sub distributeWeightedCoverage {
    my $refTree = shift;
    my $readTree = shift;

# Distribute node weights

#print "DEBUG in HLATree.pm distributeWeightedCoverage: calling _distributeNodeWeight with the the ref node ". $refTree->root->lineage. " and the Read node ". $readTree->root->lineage."\n";
    _distributeNodeWeight($refTree->root, $readTree->root);

# Recurse through the readTree and gather all the alignments and their
# associated weights
    my @alignments = @{$readTree->_returnAlignments($readTree->root)};

# when you find one, add it into the refTree
    foreach my $alignment (@alignments) {
        my ($templateLength, $imgtId, $hlaPtr) = _extractHLAfromDesc($alignment->rname);
        $refTree->_addCoverage($refTree->root,$alignment,$hlaPtr);
    }
}

## _distributeNodeWeight
# Given two Nodes, one from a reference tree and another from the read tree
# copy the weight from the read tree to the ref tree
sub _distributeNodeWeight {
    my $refNode = shift;
    my $readNode = shift;

# push the read weight from the read node onto the readWeights array of the reference node
    my @refNodeReadWeights = @{$refNode->readWeightArray};
#print "DEBUG in HLATree.pm _distributeNodeWeight: Added read weight of ".$readNode->weight." to ref readWeights array: ".join(" ", @refNodeReadWeights)."\n";
    push (@refNodeReadWeights, $readNode->weight);
    $refNode->readWeightArray(\@refNodeReadWeights);
#print "DEBUG in HLATree.pm _distributeNodeWeight: Added read weight of ".$readNode->weight." to ref readWeights array: ".join(" ", @refNodeReadWeights)."\n";

# Get children of readNode
    my @readChildren = @{$readNode->children};

# For each child of readNode, get the refNode with the same ID
    foreach my $readChild (@readChildren) {
#    print "DEBUG in HLATree.pm _distributeNodeWeight: Looking for ref node with ID ". $readChild->id." under refnode ".$refNode->lineage."\n";
        my $refChild = $refNode->getChild($readChild->id);
#    print "DEBUG in HLATree.pm _distributeNodeWeight: Found a refChild with lineage".$refChild->lineage." which matches read child with lineage ".$readChild->lineage."\n";
        _distributeNodeWeight($refChild,$readChild);

    }

# Set weight of the node
# Call _distributeNodeWeight on the located children

}

## _returnAlignments
# Given a node, returns all the alignments under that node
sub  _returnAlignments{
    my $self = shift;
    my $currentNode = shift;

    my @children = @{$currentNode->children};
    my @alignments = @{$currentNode->alignments};

    foreach my $child (@children) {
        my ($childAlignmentPtr) = $self->_returnAlignments($child);
        push (@alignments, @$childAlignmentPtr) if (scalar @$childAlignmentPtr);
    }

    return (\@alignments);
}

## _addCoverage
# Given a node and an alignment, add the weighted coverage to the correct progeeny node
# by recursing down the tree
sub _addCoverage {
    my $self = shift;
    my $currentNode = shift;
    my $alignment = shift;
    my $hlaPtr = shift;

# First check if there are any more ID's in the HLAId($hlaPtr) array
    if (scalar(@$hlaPtr)) {
        my $childId = shift @$hlaPtr;
## Recurse through children nodes until you get to a leaf
# check if any child exists with the same ID
# if it does, call _addCoverage  with a shifted HLAID array on the child
        if (my $child = $currentNode->getChild($childId)) {
            $self->_addCoverage($child, $alignment, $hlaPtr);
        }
        else {
            die "Error in HLATree.pm _addCoverage: Could not find the child node for alignment: ".$alignment->toString;
        }
    }
# Else, we're at a leaf and we can add the coverage defined in the alignment to the node
    else {
# First generate the coverage array to be added 
        my @coverage = @{$alignment->getCoverageArray};
# $currentNode->addCoverage(seqName, coveragearray)
        $currentNode->addCoverage($alignment->rname,\@coverage);
        $currentNode->addReadWeight($alignment->rname,$alignment->weight);
    }
}


#### This was a simple test method. It should be deprecated
#sub buildTreeFromFasta {
#    my $self = shift;
#    my $file = shift;
#    my $root;
#    my $seqIO = Bio::SeqIO->new(-file=> $file, -format => "fasta");
#
#    while (my $seqObj = $seqIO->next_seq()) {
#        $root = $self->root;
#        $root->_addNode($root,
#                        _extractHLAfromSimpleDesc($seqObj->display_id()),
#                        $seqObj->seq());
#    }
#}

### This was a test method using test data that no longer exists. Should be removed once
## the entire thing is working
#sub _extractHLAfromSimpleDesc {
##    my $self = shift;
#    my $desc = shift;
#    my @hla = split ":", $desc;
#    return \@hla;
#}


### This is deprecated. Should be replaced by _extractHLAfromDesc
#sub _extractHLAfromSam {
#    my $alignment = shift;
#    my @hla;
#    my @split = split "\t", $alignment;
#
#    my $refName = $split[2];
#
#    my ($imgtAcc, $HLA_type, $refLength) = $refName =~ /HLA:(HLA\d*)_(.*)_(\d+)_bp/;
#    @hla = split /:/, $HLA_type;
#    return \@hla
#}

#### The following subroutines are for testing mostly. These can be deleted eventually

# depthFirstSearchTest
sub depthFirstSearchTest {
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};

    print $currentNode->id."\n";
    foreach my $child (@children) {
        $self->depthFirstSearchTest($child);
    }
}

# depthFirstSearchTest
sub depthFirstSearchTestSMMQ {
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};

    print $currentNode->lineage."\t";
    print $currentNode->smmq."\n";
    foreach my $child (@children) {
        $self->depthFirstSearchTestSMMQ($child);
    }
}

# printWeights
sub printWeights{
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};

    print $currentNode->lineage."\t";
    print $currentNode->weight."\n";
    foreach my $child (@children) {
        $self->printWeights($child);
    }
}


sub printAlignments{
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};
    my @alignments = @{$currentNode->alignments};

    if (scalar @alignments) {
        foreach my $alignment(@alignments) {
            print $currentNode->lineage."\n";
            $alignment->print;
        }
    }

    foreach my $child (@children) {
        $self->printAlignments($child);
    }
}



sub printCoverageSums {
    my $self = shift;
    my $currentNode = shift;
    my $sum;

    my @children = @{$currentNode->children};
    my %coverages= %{$currentNode->coverage};

    if (keys %coverages) {
        foreach my $covName (keys %coverages) {
            $sum += sum( @{$coverages{$covName}});
        }
    }
    foreach my $child (@children) {
        $sum += $self->printCoverageSums($child);
    }
    print $sum."\t".$currentNode->lineage."\n";
    return $sum
}

sub printLeafCoverageSums{
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};
    my %coverages= %{$currentNode->coverage};

    if (keys %coverages) {
        foreach my $covName (keys %coverages) {
            print sum( @{$coverages{$covName}});
            print "\t$covName ".$currentNode->lineage.":::".$currentNode->id."\n";
        }
    }

    foreach my $child (@children) {
        $self->printLeafCoverageSums($child);
    }
}



sub printCoverages{
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};
    my %coverages= %{$currentNode->coverage};

    if (keys %coverages) {
        foreach my $covName (keys %coverages) {
            print "$covName ";
            print join " ", @{$coverages{$covName}};
            print "\n";
        }
    }

    foreach my $child (@children) {
        $self->printCoverages($child);
    }
}

sub printNodeWeights {
    my $self = shift;
    my $currentNode = shift;
    my %readWeights= %{$currentNode->readWeights};
    my @children = @{$currentNode->children};

    print $currentNode->lineage.":$covName\t";
    print $currentNode->weight . "\n";

    foreach my $child (@children) {
        $self->printNodeWeights($child);
    }
}

sub printReadWeightArray{
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};
    my %readWeights= %{$currentNode->readWeights};
    my @currentNodeReadWeightArray = @{$currentNode->readWeightArray};

    print join " ", ($currentNode->lineage), @currentNodeReadWeightArray;
    print "\n";

    if (keys %readWeights) {
        foreach my $covName (keys %readWeights) {
            print $currentNode->lineage.":$covName ";
            print join " ", @{$readWeights{$covName}};
            print "\n";
        }
    }

    foreach my $child (@children) {
        $self->printReadWeights($child);
    }
}




sub printSeqs{
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};
    my %sequences = %{$currentNode->seq};

    if (keys %sequences) {
        foreach my $seqName (keys %sequences) {
            print ">$seqName ".$currentNode->lineage.":::".$currentNode->id."\n".$sequences{$seqName}."\n";
        }
    }

    foreach my $child (@children) {
        $self->printSeqs($child);
    }
}


# Checks to see if any nodes with children also have sequences
sub checkNodeSeqThing {
    my $self = shift;
    my $currentNode = shift;
    my @children = @{$currentNode->children};

    if (scalar(@children) && keys %{$currentNode->{SEQ}}) {
        print "Found a node, ".$currentNode->lineage.", that has children and sequences\n";
    }
    foreach my $child (@children) {
        $self->checkNodeSeqThing($child);
    }
}

# Given an integer tier level, return all the nodes of that tier by recursion
sub get_nodes_in_tier{
    my $self  = shift;
    my $currentNode = shift;
    my $tier = shift;
    my @nodes;


    if ($tier == 0) {
        push @nodes, $currentNode;
#        print STDERR "Adding " . $currentNode->lineage . " to nodes_in_tier\n";
#        print STDERR "There are now ". scalar(@nodes) . " nodes to be returned\n";
    }
    elsif($tier > 0) {
        $tier--;
        my @children = @{$currentNode->children};
        foreach my $child (@children) {
            my @returnable_nodes = @{$self->get_nodes_in_tier($child, $tier)};
#            print STDERR "got " . scalar(@returnable_nodes) ." from function call\n";
            push @nodes, @returnable_nodes;
#            print STDERR "returning " . scalar(@nodes) ." nodes\n";
        }
    }
    else {
        die("Tier should not go negative\n");
    }

#    print STDERR "returning " . scalar(@nodes) . " \n";
    return \@nodes;
}

# Updates a crossAlignHash with the read weight between two reference nodes
sub tier_weights {
    my $self = shift;
    my $crossAlignHashPtr = shift;
    my $tier = shift;
    # first get all the nodes at the tiers
#    print STDERR "tier_weights called\n";
    my @nodes = @{$self->get_nodes_in_tier($self->root, $tier);};
    my %weights;
#    print STDERR scalar(@nodes) . " returned \n";

    # get the nodes weights
    foreach my $node (@nodes) {
        $weights{$node->lineage} = $node->weight;
    }

    print STDERR join ("\t", keys(%weights)) ." returned\n";

    my @ref_names = keys %weights;
    foreach my $ref (@ref_names) {
        foreach my $ref_again(@ref_names) {
            # Return the product of the weights of the references
            #$crossAlignHashPtr->{$ref}->{$ref_again}+=$weight{$ref_again} * $weight{$ref};
            # Return just the second reference's weight
            $crossAlignHashPtr->{$ref}->{$ref_again}+=$weights{$ref};

#        print STDERR "adding weight " .$weight{$ref_again}. " to $ref, $ref_again\n";
        }
    }
    return $crossAlignHashPtr;
}


# prune
# Given a node name, prune off the branch that contains the node name
# and recalculate read weights. This function was modified from _insertRefNode
sub prune{
    my $self = shift;
    my $currentNode = shift;
    my $hla_id_split = shift;
#    print STDERR "Prune called on ".$currentNode->lineage. " with ";
#    print join ":", @$hla_id_split, "\n";

    if (scalar(@$hla_id_split)) {
        my $childId = shift @$hla_id_split;
#        print STDERR "Looking for $childId\n";

        # Check if the child exists in this alignmentSet
        if (my $child = $currentNode->getChild($childId)) {
            $currentNode->removeAllButChild($childId);
            # Check if we should recurse further down the tree
            if(scalar(@$hla_id_split)) {
#            print STDERR "Found $childId, continuing traversal\n";
                $self->prune($child, $hla_id_split);
            }
        }
    }
}

# prune_keep_only_given
# Given an array of node nodes, keep only nodes that are specified in the list
sub prune_keep_only_given {
    my $self = shift;
    my $hla_ids_ptr = shift;
    my $debug = shift;
    my @hla_ids_split_array;
    my $max_tier = 0;

    # Parse the node names we want to keep
    foreach my $hla_id (@$hla_ids_ptr) {
        my @hla_id_split = split /:/, $hla_id;
#        shift @hla_id_split;
        push @hla_ids_split_array, \@hla_id_split;
        # get the deepest tier of the HLA_ids
        $max_tier = scalar(@hla_id_split) if (scalar(@hla_id_split) > $max_tier);
    }

    for (my $tier = 1; $tier < $max_tier; $tier++) {
        # Get all nodes in this tier
        my $nodes_in_tier_ptr = $self->get_nodes_in_tier($self->root, $tier);

        my %keeper_hash;
        my %deletors;

        print "There are " . scalar(@$nodes_in_tier_ptr) ." one nodes in tier $tier\n" if $debug;
        foreach my $node (@$nodes_in_tier_ptr) {
            my $node_lineage = $node->lineage;
            foreach my $keeper (@$hla_ids_ptr) {
                #my $keeper_lineage = $keeper->lineage;
                my $keeper_lineage = $keeper;
        # Check if the node is in the list of specified nodes to keep
        # do this by looping through and doing a regex on lineage
                if ($keeper =~ /\Q$node_lineage/) {
                    # Here if an alignment is found to a top node, prune all non keeper edges
                    $keeper_hash{$node->lineage} = $node;
#                    @deletors = ();
                    print "Found a keeper! \t tree:" if $debug;
                    print $node->lineage. " keeper:". $keeper . "\n" if $debug;
                    last;
                }
                else {
                    #push @deletors, $node;
                    $deletors{$node->lineage} = $node;
                    print "Just another bum\t tree:" if $debug;
                    print $node->lineage. " keeper:". $keeper . "\n" if $debug;
                }
            }
#            unless ($keep_node) {
#                print "trying to remove" . $node->lineage." " . $node->id."\n";
#                $node->parent()->removeChild($node->id);
#            }
        }
        # Only prune the tree if there's one of its nodes is a keeper
        if (keys %keeper_hash) {
#            foreach my $node(@deletors) {
            while (my ($lineage, $node) = each %deletors) {
                unless (exists $keeper_hash{$lineage}) {
                    print "trying to remove" . $lineage." " . $node->id."\n" if $debug;
                    $node->parent()->removeChild($node->id);
                }
            }
        }
    }
}



# multi_prune
# Given an array of node names, prune off the branches that don't contains 
# one of the node names and recalculate read weights. This prunes reference
# branches based on the order of references in an array.
# ie, the first element of the array will preferentially own all the nodes
#
# Notes:
#  - This function was modified from _insertRefNode
#  - This assigns reads to the first reference in the array that is detected

sub multi_prune {
    my $self = shift;
    my $hla_ids_ptr = shift;
    foreach my $hla_id (@$hla_ids_ptr) {
        my @hla_id_split = split /:/, $hla_id;
        shift @hla_id_split;

        if ($self->node_exists($self->root, \@hla_id_split)) {
            @hla_id_split = split /:/, $hla_id;
            shift @hla_id_split;
            $self->prune($self->root, \@hla_id_split);
            $self->_calculateWeightFromSMMQ($self->root, 1);
            last;
        }
    }
}

# multi_prune_2
# Given an array of node names, prune off the branches that don't contains 
# one of the node names and recalculate read weights.
# If multiple selected references are found in a read tree, then reads are
#  distributed to a random reference (though this could be modified to
#  distribute based on some other metric like maximum node weight in a tree
#  or randomly with observed node weight as a prior probability)
# ie, the first element of the array will preferentially own all the nodes
#
# Notes:
#  - This function was modified from _insertRefNode
#  - Assigns a read tree to a random selected reference

sub multi_prune_2 {
    my $self = shift;
    my $gene_families_ptr = shift;

#    my @families_in_tree;
    my %weights; # stores the weights of top nodes indexed by hlaid
    my @max_hla_ids;
    my $max = 0;

    # This basically determines which level to prune and stores the top nodes determined by the sum of weights for a reference
    foreach my $hla_id (keys %$gene_families_ptr) {
        my @hla_id_split = split /:/, $hla_id;
        shift @hla_id_split;

        if ($self->node_exists($self->root, \@hla_id_split)) {
#            push @families_in_tree, $hla_id;
            $weights{$hla_id} = $self->get_weight_by_hla_id($self->root, \@hla_id_split);
        }
    }

    # Find the weights of highest scoring nodes
#    if (scalar @families_in_tree) {
    if (scalar keys %weights) {
        # Sort nodes by their weight in descending value
        foreach my $hla_id (sort {$weights{$b} <=> $weights{$a}} keys %weights) {
            if ($weights{$hla_id} > $max) {
                $max = $weights{$hla_id};
                push @max_hla_ids, $hla_id;
            }
            elsif ($weights{$hla_id} == $max) {
                push @max_hla_ids, $hla_id;
            }
            else {
                last;
            }
        }
    }

    # If there are multiple nodes with the same weight, pick one randomly
    if (scalar @max_hla_ids) {
#        my $keeper = $families_in_tree[int(rand(scalar @families_in_tree))];
        my $keeper = $max_hla_ids[int(rand(scalar @max_hla_ids))];
        my @hla_id_split = split /:/, $keeper;
        shift @hla_id_split;
        $self->prune($self->root, \@hla_id_split);
        $self->_calculateWeightFromSMMQ($self->root, 1);
        return (1);
    }
    else {
#        print STDERR "Warning, none of chosen references found in this read tree. Wiping out this tree's descendants\n";
#        my @vasectomy = [];
#        $self->root->children(\@vasectomy);
#        $self = new HLATree();
        return (0);
    }
}


# Given a split array of hla id's, check if the node
# exists in the hla tree
sub node_exists {
    my $self = shift;
    my $current_node = shift;
    my $hla_id_ptr = shift;

#    print STDERR "looking for @$hla_id_ptr in " . $current_node->lineage."\n";
    if (scalar @$hla_id_ptr) {
        $childId = shift @$hla_id_ptr;
        if (my $child = $current_node->getChild($childId)) {
            return($self->node_exists($child, $hla_id_ptr));
        }
        else {
            return 0;
        }
    }
    else {
        return 1;
    }
}

# Given a split array of an hla_id, returns the weight of the node
sub get_weight_by_hla_id {
    my $self = shift;
    my $current_node = shift;
    my $hla_id_ptr = shift;

#    print STDERR "looking for @$hla_id_ptr in " . $current_node->lineage."\n";
    if (scalar @$hla_id_ptr) {
        $childId = shift @$hla_id_ptr;
        if (my $child = $current_node->getChild($childId)) {
            return($self->node_exists($child, $hla_id_ptr));
        }
        else {
            return 0;
        }
    }
    else {
        return $current_node->weight;
    }
}



return(1);
