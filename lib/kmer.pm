package kmer;
use Moose;

=head
Calculate the most frequent k-word strings in an input string

DJP-S Nov-2013
=cut


has probe_seq => (
    isa => 'Str',
    is => 'rw',
);

has len_probe_seq => (
    isa => 'Int',
    is => 'rw',
);

has mismatch_limit => (
    isa => 'Int',
    is => 'rw',
);

our @probe_array = ();

sub kmer {
    my $self = shift;

    my $string = shift;
    my $k_value = shift;

    die "Usage: kmer.pl sequence-string word-length\n" unless $string && $k_value;
 
    my $start = 0;
    my $end = length( $string ) - $k_value + 1;

    my %kwords;

    for ( my $offset = $start; $offset < $end; $offset++) {
        $kwords{ substr( $string, $offset, $k_value) } ++;
        }

    my @sorted_keys = sort { $kwords{$b} <=> $kwords{$a} } keys %kwords;
    my @vals = @kwords{@sorted_keys};

    my $max_freq =  $kwords{ $sorted_keys[0] }; # most frequent k_word score;
#    die "No subsequences exist with frequency greater than 1\n" if $max_freq == 1;

    my @frequent_strs;
    foreach my $k_key ( @sorted_keys ) { # now any others at the same frequency
        if ( $kwords{ $k_key } == $max_freq ) {
            push @frequent_strs, $k_key;
        }
        else {
            last;
        }
    }
    return \@frequent_strs;
}

=head2
Implements: search_seq_mismatch
Requires:
    self
    probe sequence (string)
    mismatch_tolerance (integer)

=cut

sub search_seq_mismatch {
    my $self = shift;


    my $probe = shift;
    my $mismatch_tolerance = shift;
    my $sequence = shift;

    $self->probe_seq( $probe );
    $self->array_probe_seq(); # sets len_probe_seq
    $self->mismatch_limit( $mismatch_tolerance );

    #generate all the subsequences the same length as the probe
    my $frequent_strs = $self->kmer( $sequence, $self->len_probe_seq );

    my @sorted_strs = sort { $self->mismatches } @{$frequent_strs};

    return \@sorted_strs;
}

sub mismatches {
    my $self = shift;

    my $ham_distance_a = $self->hamming( $a );
    my $ham_distance_b = $self->hamming( $b );

    # order the list is ascending numeric order
    return $ham_distance_a - $ham_distance_b;
}

sub hamming { # input strings are assumed to be of the same length
    my $self = shift;

    my $candidate = shift;

    my $hamming = 0;
    my @clist = split //, $candidate;

    foreach my $ind ( 0..$self->len_probe_seq-1 ) {
        if ( $probe_array[$ind] ne $clist[$ind] ){
            $hamming++;
        }
    }
    return $hamming;
}

sub array_probe_seq {
    my $self = shift;
    # for efficiency, we will always access the package variable
    # rather than the attribute.
    $self->len_probe_seq( length( $self->probe_seq ) );
    @probe_array = split //, $self->probe_seq;

    return \@probe_array;
}

sub prettify_match {
    my $self = shift;
    my $candidate = shift;

    my @new_list;
    my @clist = split //, $candidate;

    foreach my $ind ( 0..$self->len_probe_seq-1 ) {
        if ( $probe_array[$ind] ne $clist[$ind] ) {
            $new_list[$ind] = lc( $clist[$ind] );
        }
        else {
            $new_list[$ind] = uc( $clist[$ind] );
        }
    }
    return join '', @new_list;
}

__PACKAGE__->meta->make_immutable;
