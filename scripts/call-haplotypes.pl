#!/usr/bin/env perl
#
## call-haplotypes.pl
# Given output of forestUtils3.pl for tiers 2, 3, 4, and 5
# output haplotypes.

use Getopt::Long;
use warnings;
use strict;

my ($tier2_file, $tier3_file, $tier4_file, $tier5_file, $answers);
my $threshold = 0.15;

GetOptions("2=s"=>\$tier2_file, 
    "3=s"=>\$tier3_file, 
    "4=s"=>\$tier4_file, 
    "5=s"=>\$tier5_file, 
    "t|threshold=f"=>\$threshold,
);

# First read tier 2
my %tier2_data = %{read_tier_file($tier2_file)};
#print "## Tier 2##\n";
# I expect there to be 1-2 gene families for each gene
my %genes = %{get_genes_from_tier(\%tier2_data)};

#foreach my $gene (keys %genes) {
#    foreach my $node (keys %{$genes{$gene}}) {
#        print "$gene\t$node\t".$genes{$gene}{$node}."\n";
#    }
#}
#print "\n";

# Look at the evidence for both calls for a gene. Expect each call to exceed a specific threshold

# Now read tier 3
my %tier3_data = %{read_tier_file($tier3_file)};
#print "## Tier 3##\n";
## I expect there to be 1-2 gene families for each gene
#%genes = %{get_genes_from_tier(\%tier3_data)};
#foreach my $gene (keys %genes) {
#    foreach my $node (keys %{$genes{$gene}}) {
#        print "$gene\t$node\t".$genes{$gene}{$node}."\n";
#    }
#}
#print "\n";


# Here, we can get a better feel for the tier 2 calls
# The max tier3 call for the minor allele should be at least near the min call of the major allele

# Now read tier 4
my %tier4_data = %{read_tier_file($tier4_file)};

#print "## Tier 4##\n";
## I expect there to be 1-2 gene families for each gene
#%genes = %{get_genes_from_tier(\%tier4_data)};
#foreach my $gene (keys %genes) {
#    foreach my $node (keys %{$genes{$gene}}) {
#        print "$gene\t$node\t".$genes{$gene}{$node}."\n";
#    }
#}
#print "\n";

# Now read tier 5
my %tier5_data = %{read_tier_file($tier5_file)};
#print "## Tier 5##\n";
## I expect there to be 1-2 gene families for each gene
#%genes = %{get_genes_from_tier(\%tier5_data)};
#foreach my $gene (keys %genes) {
#    foreach my $node (keys %{$genes{$gene}}) {
#        print "$gene\t$node\t".$genes{$gene}{$node}."\n";
#    }
#}
#print "\n";


# Call haplotypes based on previous calling
# Tier 2 is slightly special because no previous haplotypes have been called

my ($major_t2_ptr, $minor_t2_ptr) = call_genes_from_tier2(\%tier2_data);

#print "##Tier 2 calls##\n";
#foreach my $gene (sort keys %$major_ptr) {
#    print "$gene:";
#    print "\t".$major_ptr->{$gene}->{node}."\t". $major_ptr->{$gene}->{weight} if (exists $major_ptr->{$gene});
#    print "\n";
#    print "\t".$minor_ptr->{$gene}->{node}."\t". $minor_ptr->{$gene}->{weight}  if (exists $minor_ptr->{$gene});
#    print "\n";
##}

my ($major_t3_ptr, $minor_t3_ptr) = call_haplotypes_from_tier_v3($major_t2_ptr, $minor_t2_ptr, \%tier3_data);
my ($major_t4_ptr, $minor_t4_ptr) = call_haplotypes_from_tier_v3($major_t3_ptr, $minor_t3_ptr , \%tier4_data);
my ($major_t5_ptr, $minor_t5_ptr) = call_haplotypes_from_tier_v3($major_t4_ptr, $minor_t4_ptr , \%tier5_data);

# Now fill out a haplotype table

# start backwards from tier 5 and populate major and minor haplotypes
my $major_haplotypes_ptr;
my $minor_haplotypes_ptr;
my %major_haplotypes;
my %minor_haplotypes;

#foreach my $gene (keys %$major_t5_ptr) {
#    $major_haplotypes{$gene}[4] = $major_t5_ptr->{$gene}->{node};
#}
#foreach my $gene (keys %$major_t4_ptr) {
#    $major_haplotypes{$gene}[3]= $major_t4_ptr->{$gene}->{node};
#}
#foreach my $gene (keys %$major_t3_ptr) {
#    $major_haplotypes{$gene}[2]= $major_t3_ptr->{$gene}->{node};
#}
#foreach my $gene (keys %$major_t2_ptr) {
#    $major_haplotypes{$gene}[1]= $major_t2_ptr->{$gene}->{node};
#    $major_haplotypes{$gene}[0]= $gene;
#}
#
#foreach my $gene (sort keys %major_haplotypes) {
#    print "$gene: ";
#    print join "\t", @{$major_haplotypes{$gene}};
#    print "\n";
#}

$major_haplotypes_ptr = add_tier_to_haplotype($major_haplotypes_ptr, $major_t5_ptr, 3);
$major_haplotypes_ptr = add_tier_to_haplotype($major_haplotypes_ptr, $major_t4_ptr, 2);
$major_haplotypes_ptr = add_tier_to_haplotype($major_haplotypes_ptr, $major_t3_ptr, 1);
$major_haplotypes_ptr = add_tier_to_haplotype($major_haplotypes_ptr, $major_t2_ptr, 0);

$minor_haplotypes_ptr = add_tier_to_haplotype($minor_haplotypes_ptr, $minor_t5_ptr, 3);
$minor_haplotypes_ptr = add_tier_to_haplotype($minor_haplotypes_ptr, $minor_t4_ptr, 2);
$minor_haplotypes_ptr = add_tier_to_haplotype($minor_haplotypes_ptr, $minor_t3_ptr, 1);
$minor_haplotypes_ptr = add_tier_to_haplotype($minor_haplotypes_ptr, $minor_t2_ptr, 0);


foreach my $gene (sort keys %$major_haplotypes_ptr) {
    print "ROOT:$gene\t";
    print join "\t", @{$major_haplotypes_ptr->{$gene}};
    print "\n";
    print "ROOT:$gene\t";
    if (scalar @{$minor_haplotypes_ptr->{$gene}}) {
        print join "\t", @{$minor_haplotypes_ptr->{$gene}} ;
        print "\n";
    }
    else {
        print join "\t", @{$major_haplotypes_ptr->{$gene}} ;
        print "\n";
    }
}







exit(1);


sub add_tier_to_haplotype{
    my $haplotype_ptr = shift;
    my $tier_ptr = shift;
    my $tier = shift;

    foreach my $gene (keys %$tier_ptr) {
        #$haplotype_ptr->{$gene}->[$tier]= $tier_ptr->{$gene}->{node};
        if (exists $tier_ptr->{$gene}->{node}) {
            unshift @{$haplotype_ptr->{$gene}}, $tier_ptr->{$gene}->{node};
        }
        elsif (exists  $haplotype_ptr->{$gene}->[0]){
            unshift @{$haplotype_ptr->{$gene}}, get_parent($haplotype_ptr->{$gene}->[0]);

        }
    }
    return $haplotype_ptr;
}

sub read_tier_file {
    my $tier_file = shift;
    my %data;

    open(my $tier_fh, "<", $tier_file);

    while (<$tier_fh>) {
        chomp ($_);
        my ($node_name, $value) = split /\t/, $_;
        $data{$node_name} = $value;
    }
    return \%data;
}

sub get_genes_from_tier {
    my $tier_data_ptr = shift;
    my %tier_data = %{$tier_data_ptr};
    my @references = keys %tier_data;
    my %genes;
    foreach my $reference(@references) {
        my $gene = gene_from_hla_id($reference);
        $genes{$gene}{$reference} = $tier_data{$reference};
    }
    return \%genes;
}

sub gene_from_hla_id {
    my $hla_id = shift;
    my @hla_array = split /:/, $hla_id;

    return $hla_array[1];
}

sub call_genes_from_tier2{
    my $tier_data_ptr = shift;
    my %tier_data = %{$tier_data_ptr};

    my %major_haplotypes;
    my %minor_haplotypes;

    my @nodes_sorted_by_weight = sort {$tier_data{$b} <=> $tier_data{$a}} keys %tier_data;

    foreach my $node (@nodes_sorted_by_weight) {
        my $gene = gene_from_hla_id($node);
        if (exists $major_haplotypes{$gene}) {
            if ($tier_data{$node} / $major_haplotypes{$gene}{weight} > $threshold) {
                unless (exists $minor_haplotypes{$gene}) {
                    $minor_haplotypes{$gene}{weight} = $tier_data{$node};
                    $minor_haplotypes{$gene}{node} = $node;
                }

            }
            else {
                print STDERR  "$node with ".$tier_data{$node}." did not pass the threshold\n";

            }
        }
        else {
            $major_haplotypes{$gene}{weight} = $tier_data{$node};
            $major_haplotypes{$gene}{node} = $node;
        }
    }
    return (\%major_haplotypes, \%minor_haplotypes);
}

sub call_haplotypes_from_tier_v3 {
    my $major_ptr = shift;
    my $minor_ptr = shift;
    my $data_ptr = shift;

    # foreach of the previous major and minor haplotypes, build a hash which contains the names and values of all descendants of previously called haplotypes
    my %major_descendants;
    my %minor_descendants;
    my %orphan_descendants;

    my %major_haplotypes;
    my %minor_haplotypes;

    # Basically store all the 
    foreach my $node (sort keys %$data_ptr) {
        my $gene = gene_from_hla_id($node);
        my $parent = get_parent($node);
        if (exists $major_ptr->{$gene}->{node} && $major_ptr->{$gene}->{node} eq $parent) {
#            print "Found major $node under $parent for gene $gene\n";
            $major_descendants{$gene}{$node}{node} = $node;
            $major_descendants{$gene}{$node}{weight} = $data_ptr->{$node};
        }
        elsif (exists $minor_ptr->{$gene}->{node} && $minor_ptr->{$gene}->{node} eq $parent) {
#            print "Find minor $node under $parent for gene $gene\n";
            $minor_descendants{$gene}{$node}{node} = $node;
            $minor_descendants{$gene}{$node}{weight} = $data_ptr->{$node};
        }
        else {
#            print "Could not find $node under $parent for gene $gene\n";
            $orphan_descendants{$gene}{$node}{node} = $node;
            $orphan_descendants{$gene}{$node}{weight} = $data_ptr->{$node};
        }
    } 

    # Get the largest major, minor and no descendants

    foreach my $gene(keys %minor_descendants) {
        my @sorted_minor_nodes = sort {$minor_descendants{$gene}{$b}{weight} <=> $minor_descendants{$gene}{$a}{weight}} keys %{$minor_descendants{$gene}};

        my $top_minor_node = $sorted_minor_nodes[0];
#        print "Minor" .$top_minor_node . "\t" . $minor_descendants{$gene}{$top_minor_node}{weight}."\n";
        $minor_haplotypes{$gene}{node} = $top_minor_node;
        $minor_haplotypes{$gene}{weight} = $minor_descendants{$gene}{$top_minor_node}{weight};
    }

    foreach my $gene(keys %major_descendants) {
#        print "$gene\n";
        my @sorted_major_nodes = sort {$major_descendants{$gene}{$b}{weight} <=> $major_descendants{$gene}{$a}{weight}} keys %{$major_descendants{$gene}};
        my $top_major_node = $sorted_major_nodes[0];
#        print "Major" . $top_major_node . "\t" . $major_descendants{$gene}{$top_major_node}{weight}."\n";
        $major_haplotypes{$gene}{node} = $top_major_node;
        $major_haplotypes{$gene}{weight} = $major_descendants{$gene}{$top_major_node}{weight};

        unless (exists $minor_haplotypes{$gene}) {
            if (exists $minor_ptr->{$gene}->{node}) {
                next;

            }
            if (exists $sorted_major_nodes[1]) {
                my $second_best_major_node = $sorted_major_nodes[1];
#                print "Setting second best major node as minor node $second_best_major_node\n";
                if ($major_descendants{$gene}{$second_best_major_node}{weight} > $threshold * $major_haplotypes{$gene}{weight} ) {
                    $minor_haplotypes{$gene}{node} = $second_best_major_node;
                    $minor_haplotypes{$gene}{weight} = $major_descendants{$gene}{$second_best_major_node}{weight};
#                    print "Setting no ancestor $second_best_major_node minor node";
                }
                else {
#                    print "No ancestral node $second_best_major_node did not exceed threshold\n";

                }
            }
        }

    }


    ### Don't think I need this ###
#    foreach my $gene(keys %orphan_descendants) {
#        print "$gene\n";
#        my @sorted_orphan_nodes = sort {$orphan_descendants{$gene}{$b}{weight} <=> $orphan_descendants{$gene}{$b}{weight}} keys %{$orphan_descendants{$gene}};
#        my $top_orphan_node = $sorted_orphan_nodes[0];
#        print "None ". $top_orphan_node . "\t" . $orphan_descendants{$gene}{$top_orphan_node}{weight}."\n";
#        # Check if there is already a minor node
#        unless (exists $minor_haplotypes{$gene}) {
#            if ($orphan_descendants{$gene}{$top_orphan_node}{weight} > $threshold * $major_haplotypes{$gene}{weight} ) {
#                $minor_haplotypes{$gene}{node} = $top_orphan_node;
#                $minor_haplotypes{$gene}{weight} = $orphan_descendants{$gene}{$top_orphan_node}{weight};
#                print "Setting no ancestor $top_orphan_node minor node";
#            }
#            else {
#                print "No ancestral node $top_orphan_node did not exceed threshold\n";
    #
#            }
#        }
#    }


    return (\%major_haplotypes, \%minor_haplotypes);

}
sub get_parent {
    my $refname = shift;
    my @ref_split = split /:/, $refname;
    pop @ref_split;

    return join ":", @ref_split;
}

sub call_haplotypes_from_tier_v2 {
    my $major_ptr = shift;
    my $minor_ptr = shift;
    my $data_ptr = shift;

    my %data = %$data_ptr;

    my %major_max;
    my %minor_max;

    # Ok so I get a single data point which is a node and its weight


    foreach my $node (keys %data) {
        my $gene = gene_from_hla_id($node);
        # Check to make sure 

        # initialize major_max
        unless (exists $major_max{$gene}) {
            $major_max{$gene}{weight}=0;
        }
        if (node_is_descendant($node, $major_ptr))  {
            # If it's a descendant of the major allele
            if ($data{$node} > $major_max{$gene}{weight}) {
                #            is there a minor haplotype?
                #            yes: do nothing
                unless (exists $minor_max{$gene}) {
                    #            no: Set the previous max as the minor haplotype
                    $minor_max{$gene}{weight} = $major_max{$gene}{weight};
                    $minor_max{$gene}{node} = $major_max{$gene}{node};
                }
                #       yes: set the max weight for major haplotype to this node
                $major_max{$gene}{weight} = $data{$node};
                $major_max{$gene}{node} = $node;
            }
            else {
                #       no: is there a minor haplotype?
                #           yes: do nothing
                #           no: set the minor haplotype to this node
                unless (exists $minor_max{$gene}) {
                    #            no: Set the previous max as the minor haplotype
                    $minor_max{$gene}{weight} = $major_max{$gene}{weight};
                    $minor_max{$gene}{node} = $major_max{$gene}{node};
                }
            }
        }
        # If it's a descendant of the minor allele
        elsif (node_is_descendant($node, $minor_ptr))  {
            #   is the weight the max for this desendant?
            if ($data{$node} > $minor_max{$gene}{weight}) {
                #       yes: set the max weight for the minor haplotype to this node
                $minor_max{$gene}{weight} = $data{$node};
                $minor_max{$gene}{node} = $node;
            }
            #       no: do nothing
        }
    }

    return (\%major_max, \%minor_max);
}

sub call_haplotypes_from_tier{
#    my ($major_ptr, $minor_ptr, $data_ptr) = shift;
    my $major_ptr = shift;
    my $minor_ptr = shift;
    my $data_ptr = shift;
    my %data = %{$data_ptr};
# Ok given our major and minor calls, call tier 3 haplotypes
# If no minor allele exists, see if we can call a minor allele from this data
    my %major_max;
    my %minor_max;
    foreach my $node (keys %data) {
        # if node is descendant of major allele
        my $gene = gene_from_hla_id($node);
        if (node_is_descendant($node, $major_ptr))  {

            # if node is max, set major allele to that allele
            if (exists $major_max{$gene}{weight}) {
                if ($data{$node} > $major_max{$gene}{weight}) {
                    unless (exists $minor_max{$gene}) {
                        # if node is is greater than previous max and no minor allele exists, set previous max node to minor allele
                        $minor_max{$gene}{weight} = $major_max{$gene}{weight};
                        $minor_max{$gene}{node} = $major_max{$gene}{node};
                    }
                    $major_max{$gene}{weight} = $data{$node};
                    $major_max{$gene}{node} = $node;
                }
                elsif  (! $major_max{$gene}{node} eq $node) {
                    unless (exists $minor_max{$gene}) {
                        # if node is is greater than previous max and no minor allele exists, set previous max node to minor allele
                        $minor_max{$gene}{weight} = $major_max{$gene}{weight};
                        $minor_max{$gene}{node} = $major_max{$gene}{node};
                    }
                }

            }
            else {
                $major_max{$gene}{weight} = $data{$node};
                $major_max{$gene}{node} = $node;
            }
        }
        elsif (node_is_descendant($node, $minor_ptr))  {
            if (exists $minor_max{$gene}{weight}) {
                if ($data{$node} > $minor_max{$gene}{weight}) {
                    $minor_max{$gene}{weight} = $data{$node};
                    $minor_max{$gene}{node} = $node;
                }
            }
            else {
                $minor_max{$gene}{weight} = $data{$node};
                $minor_max{$gene}{node} = $node;
            }
        }
    }

    return (\%major_max, \%minor_max);
}

sub node_is_descendant {
    my $node = shift;
    my $major_ptr = shift;
    my %major_haplotype = %$major_ptr;

    my @node_split = split /:/, $node;
    pop @node_split;
    my $parent_node = join ":", @node_split;
#    print "checking for parent of $node as ". $parent_node."\n";

    ### Debug ###
    #
    #
#    print "PARENT NODE NOT FOUND!!\n" unless($parent_node);
#    print "GENE FROM HLA ID $node NOT FOUND!!\n" unless(gene_from_hla_id($node));
#    print "major haplotype not found! $node ".gene_from_hla_id($node)." \n" unless($major_haplotype{gene_from_hla_id($node)}{node});
    ### End Debug ###
    return 0 unless (exists $major_haplotype{gene_from_hla_id($node)}{node});
    if ($major_haplotype{gene_from_hla_id($node)}{node} eq $parent_node) {
        return 1;
    }
    else {
        return 0;
    }
}
