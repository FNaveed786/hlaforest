# Package Node
# Just a generic node that has a pointer
# to parent and an array of pointers to
# children
package Node;

sub new {
    my $type = shift;
    my $self = {};
    $self->{PARENT}=undef;
    $self->{CHILDREN}=();

    bless($self,$type);
}

sub respond {
    my $self= shift;
    print "hello i am a Node\n";
}

sub parent {
    my $self = shift;
    if (@_) {$self->{PARENT} = shift}
    return $self->{PARENT};
}

sub addChild {
    my $self = shift;
    if (@_) {
        my $newChildNode = shift;
        my $childrenPtr = $self->{CHILDREN};
        push (@{$childrenPtr}, $newChildNode);
        my @sortedChildren = sort {$a->id <=> $b->id} @{$childrenPtr};

        $self->{CHILDREN} = \@sortedChildren;
    }
    else {
        print STDERR "ERROR: addChild called in Node.pm without an argument\n";
    }
}


sub children {
    my $self = shift;
# If there are extra objects being sent
    if (@_) {
        my $newChildrenPtr = shift;
        my $last = undef; 
#        print "DEBUG in Node.pm children: Adding @$newChildrenPtr to ". @{$self->children()}."\n";
# merge the new and old children arrays
        my @new = grep { ($last ne $_) && ($last = $_) } sort(@{$self->children()}, @$newChildrenPtr); 
    }
    return $self->{CHILDREN};
}

return(1);
