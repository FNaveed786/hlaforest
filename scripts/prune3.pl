#!/usr/bin/env perl

use warnings;
use strict;
use HLATree;
use Storable;
use Getopt::Long;


my @forest_files;
#my $tier;
my $debug=1;
my $threshold = 0.1;

GetOptions(
    "f=s"=>\@forest_files, 
#    "t|tier=i"=>\$tier,
    "threshold=f"=>\$threshold,
);

@forest_files = split /,/, join ",", @forest_files;




# prune3. This version should be smarter. Pruning should be done on the complete tree everytime. There is more information encoded in an unpruned tree than in a pruned tree.

# Initialize a set of selected haplotypes
my @selected_haplotypes;
my %haplotypes_hash;

# for tiers 2..5 (Gene level to Intron level)
foreach my $tier (2..5) {

# Read in a forest

# Prune i there are some selected haplotypes
#   This prune should prune iteratively then reweight once top down pruning is complete
    if (scalar @selected_haplotypes) {
        my @tmp_selected_haplotypes;
        my %genes;
        # Count the number of haplotypes for each gene
        foreach my $haplotype(@selected_haplotypes) {
            my $gene = gene_from_hla_id($haplotype);
            $genes{$gene}++;
        }

#        print join "\n", @selected_haplotypes, "";

        my ($pruned_weight_hash_ptr) = prune(\@forest_files, \@selected_haplotypes, $tier, 0);
        my @sorted_haplotypes = @{sort_keys_desc($pruned_weight_hash_ptr)};

#        print "Tier $tier\n";
#        print_hash_ptr($pruned_weight_hash_ptr);

#        # Call haplotypes of selected
#        foreach my $haplotype(@selected_haplotypes) {
##            print "Looking for child of $haplotype\n";
#            foreach my $selected_haplotype (@sorted_haplotypes) {
#                if($selected_haplotype =~ /\Q$haplotype/) {
#                    push @tmp_selected_haplotypes, $selected_haplotype;
#                    last;
#                }
#            }
#        }
#
#        $pruned_weight_hash_ptr = prune(\@forest_files, \@tmp_selected_haplotypes, $tier, 0);
#        print "Tier $tier.5\n";
#        print_hash_ptr($pruned_weight_hash_ptr);
#        # For all genes with only 1 haplotype, check if the pruned haplotype exceeds the threshold




        foreach my $haplotype(@selected_haplotypes) {
        # if the gene only has a single haplotype select two if it passes the threshold
            if ($genes{gene_from_hla_id($haplotype)}==1) {
#                print "Hey now, i only have one of this gene\n";
                my $num_selected_haplotypes = 0;
                my $first_haplotype;
                foreach my $selected_haplotype (@sorted_haplotypes) {
                    last if ($num_selected_haplotypes > 1);
#                    print "is there a $haplotype in $selected_haplotype?\n";
                    if($selected_haplotype =~ /\Q$haplotype/) {
                        if ($first_haplotype) {
                            if ($pruned_weight_hash_ptr->{$first_haplotype} / $pruned_weight_hash_ptr->{$selected_haplotype}) {
#                                print "found one baby $selected_haplotype\n";
                                push @tmp_selected_haplotypes, $selected_haplotype;
                                $num_selected_haplotypes++;
                            }
                            else {
                                last;
                            }
                        }
                        else {
                            $first_haplotype = $selected_haplotype;
                            push @tmp_selected_haplotypes, $selected_haplotype;
                            $num_selected_haplotypes++;
                        }
                    }
                }
            }
            # Pick the child for each haplotype
            elsif ($genes{gene_from_hla_id($haplotype)}==2) {
                foreach my $selected_haplotype (@sorted_haplotypes) {
                    if($selected_haplotype =~ /\Q$haplotype/) {
                        push @tmp_selected_haplotypes, $selected_haplotype;
                        last;
                    }
                }
            }
        else {
            print "I'm not sure what to do here. I love lamp\n";
            exit(0);
        }
    }
        @selected_haplotypes = @tmp_selected_haplotypes;
        %haplotypes_hash = %{add_haplotypes_to_hash(\%haplotypes_hash, \@selected_haplotypes)};

    }

# Select the first round
    else {
        my ($sum_hash_ptr) = get_weights(\@forest_files, $tier);
#        print_hash_ptr($sum_hash_ptr);

# Get the top node
        my ($top_nodes_hash_ptr) = get_node_weights($sum_hash_ptr, $tier);

# # Select the top child for each node at the tier (haplotype 1)
        my ($ordered_top_node_lineages_array_ptr) = top_weighted_nodes_hash_to_array ( $top_nodes_hash_ptr );

# Prune
        my ($pruned_weight_hash_ptr) = prune(\@forest_files, $ordered_top_node_lineages_array_ptr, $tier);
#        print_hash_ptr($pruned_weight_hash_ptr);

        my ($top_children_hash_ptr) = get_top_2_children($pruned_weight_hash_ptr, $tier, $threshold);
#        print_hash_ptr($top_children_hash_ptr);
        @selected_haplotypes = keys %$top_children_hash_ptr;
        %haplotypes_hash = %{add_haplotypes_to_hash(\%haplotypes_hash, \@selected_haplotypes)};
    }

# Get the second node at the tier
# Check the threshold and store the selected haplotypes 
# if they exceed the threshold
# endfor
}

hash_walk(\%haplotypes_hash, [], \&print_keys);



exit(1);


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
sub prune {
    my $forest_files_ptr = shift;
    my $ordered_lineages_ptr = shift;
    my $tier = shift;
    my $debug = shift;

    my %prune_weights;
    foreach my $forest_file ( @forest_files ) {
        my $hla_forest_ptr = eval { retrieve( $forest_file ) };

        foreach my $hla_tree (@{$hla_forest_ptr}) {
            print "New tree\n" if $debug;
            $hla_tree->prune_keep_only_given($ordered_lineages_ptr, $debug);
            $hla_tree->_calculateWeightFromSMMQ($hla_tree->root, 1);
            my @nodes = @{$hla_tree->get_nodes_in_tier($hla_tree->root, $tier)};
            foreach my $node (@nodes) {
                my $weight = $node->weight;
                print "Adding weight for ". $node->lineage."\n" if $debug;
                if ($weight) {
                    $prune_weights{$node->lineage} += $node->weight;
                }
                else {
                    print $node->id . "did not have a weight\n" if $debug;
                }
            }
        }
    }
    return \%prune_weights;
}

# Given a weight hash, returns an array containing the top two children of 
# all nodes in tier 2
sub get_top_2_children{ 
    my $weights_hash_ptr = shift;
    my $tier = shift;
    my $threshold = shift;

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








##### Helpers #####
# Helper function to print the contents of a hash pointer
sub print_hash_ptr {
    my $hash_ptr = shift;
    foreach my $key ( sort {$hash_ptr->{$b} <=> $hash_ptr->{$a}} keys %$hash_ptr) {
        print join "\t", $key, $hash_ptr->{$key}, "\n";
    }
}

# return sorted array of keys
sub sort_keys_desc {
    my $hash_ptr = shift;
    my @sorted_keys;
    foreach my $key ( sort {$hash_ptr->{$b} <=> $hash_ptr->{$a}} keys %$hash_ptr) {
        push @sorted_keys, $key;
    }
    return \@sorted_keys;
}


# Given an hla id, return the gene
sub gene_from_hla_id {
    my $hla_id = shift;
    my @hap_split = split /:/, $hla_id;
    return $hap_split[1];
}

# Taken from http://stackoverflow.com/questions/11505100/perl-how-to-turn-array-into-nested-hash-keys
sub add_haplotypes_to_hash {
    my $hash_ptr = shift;
    my $haplotypes_ptr = shift;
    foreach my $haplotype (@$haplotypes_ptr) {
        my @haplotype_split = split /:/, $haplotype;
        $hash_ptr = array_to_hash_keys($hash_ptr, @haplotype_split);
    }
    return $hash_ptr;
}


sub array_to_hash_keys {
    my $ref   = \shift;  
    my $h     = $$ref;
    $ref      = \$$ref->{ $_ } foreach @_;
    return $h;
}

# Used to print the hash
# from http://stackoverflow.com/questions/2363142/how-to-iterate-through-hash-of-hashes-in-perl
sub hash_walk {
    my ($hash, $key_list, $callback) = @_;
    while (my ($k, $v) = each %$hash) {
        # Keep track of the hierarchy of keys, in case
        # our callback needs it.
        push @$key_list, $k;

        if (ref($v) eq 'HASH') {
            # Recurse.
            hash_walk($v, $key_list, $callback);
        }
        else {
            # Otherwise, invoke our callback, passing it
            # the current key and value, along with the
            # full parentage of that key.
            $callback->($k, $v, $key_list);
        }

        pop @$key_list;
    }
}

sub print_keys{
    my ($k, $v, $key_list) = @_;
#    printf "k = %-8s  key_list = [%s]\n", $k, "@$key_list";
    print join ":", @$key_list;
    print "\n";
}
