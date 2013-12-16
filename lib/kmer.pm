package kmer;
use Moose;

=head agtc
Calculate the most frequent k-word strings in an input string.

Handles mismatches so that word-length substrings can be returned with up to
mismatch_value mismatches distributed across the length of the string. 

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

has title => (
    isa => 'Str',
    is => 'rw',
);

our @probe_array = ();

=head kmer
Given sequence-string word-length
Returns a list of all strings of length word-length in the supplied string

=cut

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

=head kmer_mismatch
Given sequence_string probe_string mismatch_value hash_reference

Returns the list of frequently occurring strings matching the probe_string
with a maximum number of mismatches given by the mismatch_value

The final hash reference is a pointer to the hash that will contain the growing
list of candidates.
=cut

sub kmer_mismatch {
    my $self = shift;

    my $string = shift;
    my $candidates = shift;

    my $k_value = $self->len_probe_seq;
    my $mismatch_value = $self->mismatch_limit;

    my $start = 0;
    my $end = length( $string ) - $k_value + 1;


    for ( my $offset = $start; $offset < $end; $offset++) {
        my $candidate = substr( $string, $offset, $k_value);
        my $hamd = $self->hamming( $candidate );
        if ( $hamd <= $mismatch_value ) {
            $candidates->{ $candidate } ++;
        }
    }
    return ;
}


=head seq_frequency
Given a hash of sequence keys

Returns a list of sorted keys based on hamming distance from the probe
=cut

sub seq_frequency {
    my $self = shift;
    my $kwords = shift; # ref to a hash
$DB::single=1;
# Sort the keys on their frequency (value) fields
    my @sorted_keys = sort {
        $kwords->{$b} <=> $kwords->{$a}
    } keys %{$kwords};
    my @vals = @{$kwords}{@sorted_keys};

    my $max_freq =  $kwords->{ $sorted_keys[0] }; # most frequent k_word score;
#    die "No subsequences exist with frequency greater than 1\n" if $max_freq == 1;

    my @frequent_strs;
    foreach my $k_key ( @sorted_keys ) { # now any others at the same frequency
        if ( $kwords->{ $k_key } == $max_freq ) {
            push @frequent_strs, $k_key;
        }
        else {
            last;
        }
    }
    return \@frequent_strs;
}

=head2 search_seq_mismatch
Implements: search_seq_mismatch
Requires:
    self
    probe sequence (string)
    mismatch_tolerance (integer)

Algorithm:
generate all the subsequences of the length of the probe sequence
generate a list sorted by hamming distance

This will calculate every mismatch using the perl sort algorithm

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
    my $frequent_strs = $self->kmer_mismatch( $sequence, $self->len_probe_seq );

    my @sorted_strs = sort { $self->mismatches } @{$frequent_strs};

    return \@sorted_strs;
}

=head search_seq_mismatch_file
Given a probe_sequence, mismatch_value and file_name

Return a list of candidate sequences within the mismatch threshold, sorted
from min to max hamming edit distance against the probe sequence.

The file should be in Fasta format. The first line will be used as the
title data for the sequence.

This method is designed to handle large files effectively (if not efficiently!)
=cut

sub search_seq_mismatch_file {
    my $self = shift;
$DB::single=1;
    my $probe = shift;
    my $mismatch_tolerance = shift;
    my $sequence_file = shift;

    open( my $seq_fh, "<", $sequence_file )
        or die "cannot open < $sequence_file: $!";

    my $title = <$seq_fh> ; 
    chomp $title;
    chop($title) if ($title =~ m/\r$/);
    $self->title( $title );

    $self->probe_seq( $probe );
    $self->array_probe_seq(); # sets len_probe_seq
    $self->mismatch_limit( $mismatch_tolerance );
    my $residual_length = $self->len_probe_seq - 1; # The amount of overlap required
    # to provide a seamless join for word_length words.


    #generate all the subsequences the same length as the probe
    my %candidates;
    my $residual = '';
    my $feedback = 0;
    my $line_number = 0;
    while (defined( my $line = <$seq_fh>) ) { # some loop to go through the file for portions of the sequence to send in to kmer_mismatch
        ++$line_number;
        if ( ++$feedback > 100000 ) {
            print $line_number . "\n";
            $feedback = 1;
        }
        chomp $line;
        chop($line) if ($line =~ m/\r$/);
        my $sequence = $residual . $line;
        $residual = substr( $sequence, -1 * $residual_length); # sets that up for next time round
        $self->kmer_mismatch( $sequence, \%candidates);
    }

    my $frequent_strs = $self->seq_frequency( \%candidates );
    my @sorted_strs = sort { $self->mismatches } @{$frequent_strs};

    return \@sorted_strs;
}

=head mismatches
Given globals $a and $b

Returns the hamming edit distance between the two strings
=cut

sub mismatches {
    my $self = shift;

    my $ham_distance_a = $self->hamming( $a );
    my $ham_distance_b = $self->hamming( $b );

    # order the list is ascending numeric order
    return $ham_distance_a - $ham_distance_b;
}

=head hamming
Given a candidate sequence
Returns the hamming edit distance with the probe sequence stored in the global probe_array
=cut

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
