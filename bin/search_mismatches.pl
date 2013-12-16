use strict;
use kmer;


my $probe = shift @ARGV;
my $mismatch_tolerance = shift @ARGV;
my $sequence_file_name = shift @ARGV;

my $k = kmer->new;

my $mismatch_strs = $k->search_seq_mismatch_file( $probe, $mismatch_tolerance, $sequence_file_name );

print $k->title . "\n\n";
print $probe . "(probe)\n\n";
if ( @$mismatch_strs ) {
    foreach my $match ( @$mismatch_strs ) {
        my $ham = $k->hamming( $match );
        last if  $ham > $mismatch_tolerance;
        printf $k->prettify_match( $match ) . ' => ' . $ham  . "\n";
    }
}
