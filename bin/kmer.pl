use strict;

=head
Calculate the most frequent k-word strings in an input string

DJP-S Nov-2013
=cut

my ($string, $k_value) = @ARGV;

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
die "No subsequences exist with frequency greater than 1\n" if $max_freq == 1;

my @frequent_strs;
foreach my $k_key ( @sorted_keys ) { # now any others at the same frequency
    if ( $kwords{ $k_key } == $max_freq ) {
        push @frequent_strs, $k_key;
    }
    else {
        last;
    }
}
my @sfrequent_strs = sort @frequent_strs;

if ( @sfrequent_strs ) {
    printf join( ' ', @sfrequent_strs ) . "\n";
}
