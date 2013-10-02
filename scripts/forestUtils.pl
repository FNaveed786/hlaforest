#!/usr/bin/env perl

use warnings;
use strict;
use HLATree;
use Storable;
use Getopt::Long;
use Data::Dumper;

my $forestFile;
my $sum_flag;
my $uniq_sum_flag;
my $circos_flag;
my $R_flag;
my $tier = 1;
my $out_forest;
my $prune;
my $gene_filter;
my $d3_flag;

GetOptions(
    "f=s"=>\$forestFile, 
    "t|tier:i"=>\$tier, 
    "circos"=>\$circos_flag, 
    "d3_graph"=>\$d3_flag,
    "sum"=>\$sum_flag, 
    "unique_sum"=>\$uniq_sum_flag, 
    "o|out:s"=>\$out_forest, 
    "p|prune:s"=>\$prune,
    "R"=>\$R_flag,
    "gene_filter:s"=>\$gene_filter,
);

#@forestFiles = split /,/, join ",", @forestFiles;

print STDERR "Reading in hla forest\t";
my $hla_forest_ptr = eval { retrieve( $forestFile ) };
print STDERR "done\n";

my %crossAlignHash;
my $num_reads = 0;
# print Data::Dumper->Dump(@hla_forest);
print STDERR "num of forests: ". scalar(@{$hla_forest_ptr}). "\n";

if ($circos_flag) {
    my @gene_filter_array;

    # This adds weights to a hash that stores cross alignment information
    foreach my $hla_tree (@{$hla_forest_ptr}) {
    if ($gene_filter) {
        @gene_filter_array = split /:/, $gene_filter;
    }
        # if the gene filter exists then check if the gene exists within in the tree
        if ($gene_filter and $hla_tree->node_exists($hla_tree->root, \@gene_filter_array )) {
            $hla_tree->tier_weights(\%crossAlignHash, $tier);
        }
        # if the tree doesn't have the child then skip it
        elsif ($gene_filter) {
            # skip this one
        }
        # if there's no filter specified, then don't bother filtering
        else {
            $hla_tree->tier_weights(\%crossAlignHash, $tier);
        }
    }

    my @refs = sort (keys %crossAlignHash);
    my $min_weight = 0;

    print join "\t",  "col", @{nodeIdArrayToCircosFriendlyNames(\@refs)};
    print "\n";
    foreach my $ref_one (@refs) {
        # need a circos friendly print
#        print $ref_one;
        print nodeIdToCircosFriendlyName($ref_one);
        foreach my $ref_two (@refs) {
            if (exists $crossAlignHash{$ref_one}{$ref_two}) {
                if ($min_weight) {
                    if($crossAlignHash{$ref_one}{$ref_two} > $min_weight) {
                        print "\t". $crossAlignHash{$ref_one}{$ref_two};
                    }
                    # If an edge does not exist, print a circos friend '-' instead of a '0'
                    else {
                        print "\t-";
                    }
                }
                else {
                    print "\t". $crossAlignHash{$ref_one}{$ref_two};
                }
            }
            else {
                print "\t-";
            }
        }
        print "\n";
    }
}

elsif($d3_flag) {
    my $tree = HLATree->new();

    foreach my $hla_tree (@{$hla_forest_ptr}) {
        $tree->addTreeWeights($hla_tree);
    }

    $tree->d3WeightJson($tree->root());
}

elsif ($R_flag) {
    foreach my $hla_tree (@{$hla_forest_ptr}) {
        $hla_tree->tier_weights(\%crossAlignHash, $tier);
    }

    my @refs = sort (keys %crossAlignHash);
    my $min_weight = 0;

    print join "\t",  "col", @refs;
    print "\n";
    foreach my $ref_one (@refs) {
        print $ref_one;
        foreach my $ref_two (@refs) {
            if (exists $crossAlignHash{$ref_one}{$ref_two}) {
                if ($min_weight) {
                    if($crossAlignHash{$ref_one}{$ref_two} > $min_weight) {
                        print "\t". $crossAlignHash{$ref_one}{$ref_two};
                    }
                    else {
                        print "\t0";
                    }
                }
                else {
                    print "\t". $crossAlignHash{$ref_one}{$ref_two};
                }
            }
            else {
                print "\t0";
            }
        }
        print "\n";
    }
}


elsif($sum_flag) {
    my %weights;
    foreach my $hla_tree (@{$hla_forest_ptr}) {
        my @nodes = @{$hla_tree->get_nodes_in_tier($hla_tree->root, $tier)};
        foreach my $node (@nodes) {
            $weights{$node->lineage} += $node->weight;
        }
    }

    foreach my $ref (keys %weights) {
        print join "\t", $ref, $weights{$ref};
        print "\n";
    }
}

elsif($uniq_sum_flag) {
    my %weights;
    foreach my $hla_tree (@{$hla_forest_ptr}) {
        my @nodes = @{$hla_tree->get_nodes_in_tier($hla_tree->root, $tier)};
        if (scalar(@nodes) == 1) {
            foreach my $node (@nodes) {
                $weights{$node->lineage} += $node->weight;
            }
        }
    }

    foreach my $ref (keys %weights) {
        print join "\t", $ref, $weights{$ref};
        print "\n";
    }
}



elsif($prune) {
    print STDERR "Splitting reference name $prune";
    foreach my $hla_tree (@{$hla_forest_ptr}) {
        my @hla_id_split = split ":", $prune;
        # Check if the hla_id_split exists in alignment set
        if ($hla_tree->node_exists($hla_tree->root, \@hla_id_split)) {
            @hla_id_split = split ":", $prune;
            # prune the tree
            $hla_tree->prune($hla_tree->root, \@hla_id_split);
            $hla_tree->_calculateWeightFromSMMQ($hla_tree->root, 1)
        }
    }

    my $result = eval { store( $hla_forest_ptr, $out_forest ) };

    if( $@ )
    { warn "Serious error from Storable: $@" }
    elsif( not defined $result )
    { warn "I/O error from Storable: $!" }
}

# Given an hla node in the style of 'ROOT:A:01' return a circos friendly name like 'A_01'
sub nodeIdToCircosFriendlyName {
    my $id = shift;
    my @id_split = split /:/, $id;
    shift @id_split;
    return join "_", @id_split;
}

# Given an array of HLA ids, return an array of circos friendly names
sub nodeIdArrayToCircosFriendlyNames {
    my $id_array_ptr = shift;
    my @circos_friendly_array;

    foreach my $id (@$id_array_ptr) {
        push @circos_friendly_array, nodeIdToCircosFriendlyName($id);
    }
    return \@circos_friendly_array;
}
