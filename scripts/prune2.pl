#!/usr/bin/env perl

use warnings;
use strict;
use HLATree;
use Storable;
use Getopt::Long;


my @forest_files;
my $tier;
my $debug=1;
my $threshold = 0.1;

GetOptions(
    "f=s"=>\@forest_files, 
    "t|tier=i"=>\$tier,
    "threshold=f"=>\$threshold,
);

@forest_files = split /,/, join ",", @forest_files;

pick_children(\@forest_files, $tier);


exit(1);

# These functions could be offloaded into an HLAForest.pm which also includes
# functions to read and write forests using Storable
sub pick_children {
    my $forest_files_ptr = shift;
    my $tier = shift;

# This subroutine should 
# Given a tier, return the unique sums and the multi-mapped sums
    my ($sum_hash_ptr) = get_weights($forest_files_ptr, $tier);

#    print_hash_ptr($sum_hash_ptr);

# Select the top child for each node at the tier (haplotype 1)
    my ($top_nodes_hash_ptr) = get_node_weights($sum_hash_ptr, $tier);

#    print_top_nodes_hash_ptr($top_nodes_hash_ptr);

   my ($ordered_top_node_lineages_array_ptr) = top_weighted_nodes_hash_to_array ( $top_nodes_hash_ptr );
    
#    print join "\t", sort @$ordered_top_node_lineages_array_ptr, "\n";

# Prune trees with the top children
    my ($pruned_weight_hash_ptr) = prune_and_get_weights($forest_files_ptr, $ordered_top_node_lineages_array_ptr, $tier);

#    print_hash_ptr($pruned_weight_hash_ptr);

# Select the next best top children (haplotype 2)
    my ($top_children_hash_ptr) = get_top_2_children($pruned_weight_hash_ptr, $tier);
#    print_hash_ptr($top_children_hash_ptr);

# Reclassifiy trees into one of four groups
# 1. Trees mapping uniquely to one of the top children (unique trees)
# 2. Trees mapping to multiple children including one of the top children (multi trees)
# 3. Trees mapping to none of the top children (discordant trees)
# Save trees into one of those three groups
# Output statistics for all children including: unique evidence, unpruned evidence, pruned evidence
    my ($repruned_weight_hash_ptr) = prune2_and_get_weights($forest_files_ptr, $top_children_hash_ptr, $tier);
    print_hash_ptr($repruned_weight_hash_ptr);
}


# pick_children2 aims to handle unique and multi-mapped trees separately
# not sure if this is the best method
sub pick_children2 {
    my $forest_files_ptr = shift;
    my $tier = shift;

# This subroutine should 
# Given a tier, return the unique sums and the multi-mapped sums
    my ($unique_sum_hash_ptr, $multi_sum_hash_ptr) = get_unique_and_multi_weights($forest_files_ptr, $tier);
    #print_hash_ptr($unique_sum_hash_ptr);
    #print_hash_ptr($multi_sum_hash_ptr);

# Select the top child for each node at the tier (haplotype 1)
    my ($unique_top_nodes_hash_ptr) = get_node_weights($unique_sum_hash_ptr, $tier);
    my ($multi_top_nodes_hash_ptr) = get_node_weights($multi_sum_hash_ptr, $tier);

    print_top_nodes_hash_ptr($unique_top_nodes_hash_ptr);
    print_top_nodes_hash_ptr($multi_top_nodes_hash_ptr);

    # Check to make sure that the unique and multi weight sums agree here

    #print_hash_ptr(check_top_nodes_hashes($unique_top_nodes_hash_ptr, $multi_top_nodes_hash_ptr));
    

    my ($ordered_unique_top_node_lineages_array_ptr) = top_weighted_nodes_hash_to_array ( $unique_top_nodes_hash_ptr );
    my ($ordered_multi_top_node_lineages_array_ptr) = top_weighted_nodes_hash_to_array ( $multi_top_nodes_hash_ptr );
    
    print join "\t", sort @$ordered_unique_top_node_lineages_array_ptr, "\n";
    print "\n";
    print join "\t", sort @$ordered_multi_top_node_lineages_array_ptr, "\n";




# Prune trees with the top children
# Recalculate unique and multi-mapped sums
# Select the next best top children (haplotype 2)
# Reclassifiy trees into one of four groups
# 1. Trees mapping uniquely to one of the top children (unique trees)
# 2. Trees mapping to multiple children including one of the top children (multi trees)
# 3. Trees mapping to none of the top children (discordant trees)
# Save trees into one of those three groups
# Output statistics for all children including: unique evidence, unpruned evidence, pruned evidence



}

# Given an array pointer of forest files, return two pointers to hashes containing
# the weights at each node at the specified tier
sub get_weights {
    my $forest_files_ptr = shift;
    my $tier = shift;
    my %weights;

    foreach my $forest_file ( @forest_files ) {
        my $hla_forest_ptr = eval { retrieve( $forest_file ) };

        foreach my $hla_tree (@{$hla_forest_ptr}) {
            my @nodes = @{$hla_tree->get_nodes_in_tier($hla_tree->root, $tier)};
            foreach my $node (@nodes) {
                $weights{$node->lineage} += $node->weight;
            }
        }
    }
    return (\%weights);
}



# Given an array pointer of forest files, return two pointers to hashes containing
# the weights at each node at the specified tier
sub get_unique_and_multi_weights {
    my $forest_files_ptr = shift;
    my $tier = shift;
    my %unique_weights;
    my %multi_weights;

    foreach my $forest_file ( @forest_files ) {
        my $hla_forest_ptr = eval { retrieve( $forest_file ) };

        foreach my $hla_tree (@{$hla_forest_ptr}) {
            my @nodes = @{$hla_tree->get_nodes_in_tier($hla_tree->root, $tier)};
            if (scalar(@nodes) == 1) {
                foreach my $node (@nodes) {
                    $unique_weights{$node->lineage} += $node->weight;
                }
            }
            else {
                foreach my $node (@nodes) {
                    $multi_weights{$node->lineage} += $node->weight;
                }
            }
        }
    }
    return (\%unique_weights, \%multi_weights);
}

# Given a hash of weights genereated by get_weights(), return to the top children of a specified tier
sub get_node_weights{
    my $weight_hash_ptr = shift;
    my $tier = shift;

    my %weights = %$weight_hash_ptr;
    my @sorted_ref_by_weight = sort {$weights{$a} <=> $weights{$b}} keys %weights;
    my %top_nodes;
    my $parent_tier = $tier -1;

    foreach my $ref (@sorted_ref_by_weight) {
        my @ref_split = split /:/, $ref;
        my $parent_lineage = join ':', @ref_split[1..$parent_tier];
        unless (exists $top_nodes{$parent_lineage}) {
            $top_nodes{$parent_lineage}{VALUE} = $weights{$ref};
            $top_nodes{$parent_lineage}{REFNAME} = $ref;
        }
        elsif ($weights{$ref} > $top_nodes{$parent_lineage}{VALUE}) {
            $top_nodes{$parent_lineage}{VALUE} = $weights{$ref};
            $top_nodes{$parent_lineage}{REFNAME} = $ref;
        }
    }

    return \%top_nodes;
}

# This method is essentially running 'keys' on the hash of weights generated by get_weights(), but because it's got a funny
# structure it needs some more code
sub top_weighted_nodes_hash_to_array {
    my $top_weighted_node_hash_ptr = shift;
    my %top_genes = %$top_weighted_node_hash_ptr;

    my @ordered_first_ref = sort {$top_genes{$a}{VALUE} <=> $top_genes{$b}{VALUE}} keys %top_genes;
    my @ordered_first_lineage;
    foreach my $ref (@ordered_first_ref) {
        push (@ordered_first_lineage, $top_genes{$ref}{REFNAME});
    }

    return \@ordered_first_lineage;
}

# Given an array of forest files, an array of the lineages with the most evidence, and 
# a tier, return a recalculated weight hash
sub prune_and_get_weights {
    my $forest_files_ptr = shift;
    my $ordered_lineages_ptr = shift;
    my $tier = shift;

    my %prune_weights;
    foreach my $forest_file ( @forest_files ) {
        my $hla_forest_ptr = eval { retrieve( $forest_file ) };

        foreach my $hla_tree (@{$hla_forest_ptr}) {
            $hla_tree->multi_prune($ordered_lineages_ptr);
            my @nodes = @{$hla_tree->get_nodes_in_tier($hla_tree->root, $tier)};
            foreach my $node (@nodes) {
                $prune_weights{$node->lineage} += $node->weight;
            }
        }
    }

   
    #print_hash_ptr(\%prune_weights) if $debug;

    return \%prune_weights;
}


# Given a weight hash, returns an array containing the top two children of 
# all nodes in tier 2
sub get_top_2_children{ 
    my $weights_hash_ptr = shift;
    my $tier = shift;

    my %weights = %$weights_hash_ptr;
    my $parent_tier = $tier - 1;
    my %gene_weights;
    my %top_gene_families;

    # for each haplotype sorted by weight
    foreach my $lineage(sort {$weights{$b}<=>$weights{$a}} keys %weights) {
        my @ref_split = split /:/, $lineage;
        # Push all the weights of each lineage onto a hash indexed by the parent tier of the haplotype
        my $parent_lineage = join ':', @ref_split[1..$parent_tier];
        push @{$gene_weights{$parent_lineage}}, $lineage;
    }

    foreach my $gene (keys %gene_weights) {
#    push @top_gene_families, @{$gene_weights{$gene}}[0..1];
        # Push the top scoring haplotype onto the list
        $top_gene_families{$gene_weights{$gene}[0]} = 1;
        # Push the second top scoring haplotype onto the list if it exceeds the threshold
        # And if it even exists
        if ($gene_weights{$gene}[1]) {
            if ($weights{$gene_weights{$gene}[1]} / $weights{$gene_weights{$gene}[0]} > $threshold) {
                $top_gene_families{$gene_weights{$gene}[1]} = 1;
            }
        }
    }

    return \%top_gene_families;
}

# 
sub prune2_and_get_weights {
    my $forest_files_ptr = shift;
    my $top_nodes_hash_ptr = shift;;
    my $tier = shift;

    my @forest_files = @$forest_files_ptr;
    my %top_nodes = %$top_nodes_hash_ptr;
    my %prune_weights_2;


    my $nonconforming_tree_count = 0;
    foreach my $forest_file ( @forest_files ) {
        my $hla_forest_ptr = eval { retrieve( $forest_file ) };
        my @pruned_trees;
        my @non_conformist_trees;

        foreach my $hla_tree (@{$hla_forest_ptr}) {
            if ($hla_tree->multi_prune_2(\%top_nodes)) { #this returns 1 if successfully pruned
                my @nodes = @{$hla_tree->get_nodes_in_tier($hla_tree->root, $tier)};
                foreach my $node (@nodes) {
                    $prune_weights_2{$node->lineage} += $node->weight;
                }
                push @pruned_trees, $hla_tree;
            }
            else {
                # set non conforming reads aside
                $nonconforming_tree_count++;
                push @non_conformist_trees, $hla_tree;
            }
        }
        eval {store (\@pruned_trees, $forest_file."pruned");};
        eval {store (\@non_conformist_trees, $forest_file."nonconforming");};
    }


    return \%prune_weights_2;
}


##### Helpers #####
# Helper function to print the contents of a hash pointer
sub print_hash_ptr {
    my $hash_ptr = shift;

    foreach my $key ( sort {$hash_ptr->{$b} <=> $hash_ptr->{$a}} keys %$hash_ptr) {
        print join "\t", $key, $hash_ptr->{$key}, "\n";
    }
}

# Helper function to print the contents of a top nodes hash pointer
sub print_top_nodes_hash_ptr {
    my $hash_ptr = shift;
    my %hash = %$hash_ptr;

    foreach my $key ( sort {$hash{$b}{VALUE} <=> $hash{$a}{VALUE}} keys %hash) {
        print join "\t", $key, $hash_ptr->{$key}->{REFNAME}, $hash_ptr->{$key}->{VALUE}, "\n";
    }
}

# Given two hash of weights generated by get_weights, returns the agreement between both
sub check_top_nodes_hashes {
    my $hash_a_ptr = shift;
    my $hash_b_ptr = shift;
    my %hash_a = %$hash_a_ptr;
    my %hash_b = %$hash_b_ptr;
    my %ref_count;

    foreach my $key ( keys %hash_a)  {
        $ref_count{$hash_a{$key}{REFNAME}}++;
    }
    foreach my $key ( keys %hash_b)  {
        $ref_count{$hash_b{$key}{REFNAME}}++;
    }
    return \%ref_count;
}


