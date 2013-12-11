use strict;
use kmer;

#my $sequence = 'AAACCGCAATTTAAAACCCGGCCGACTACTTTCGAT';
#my $probe = 'AAACCGC';
#my $mismatch_tolerance = 6;

my $sequence = shift @ARGV;
my $probe = shift @ARGV;
my $mismatch_tolerance = shift @ARGV;

my $k = kmer->new;

my $mismatch_strs = $k->search_seq_mismatch( $probe, $mismatch_tolerance, $sequence );

print $probe . "\n\n";
if ( @$mismatch_strs ) {
    foreach my $match ( @$mismatch_strs ) {
        my $ham = $k->hamming( $match );
        last if  $ham > $mismatch_tolerance;
        printf $k->prettify_match( $match ) . ' => ' . $ham  . "\n";
    }
}
