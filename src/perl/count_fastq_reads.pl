#!/usr/local/bin/perl
#fpd120202_2228 

use strict;
use warnings;

main() ;

sub main {

   my $fastq_fn = $ARGV[0] || die "ERROR: must specify FASTQ file name" ;
   if (!-e $fastq_fn) {
      die "file not found: $fastq_fn\n";
   } elsif (!-s $fastq_fn) {
      die "file is empty: $fastq_fn\n";
   }

   my $wc_out ;
   if ($fastq_fn =~ /gz$/) {
      $wc_out = `zcat $fastq_fn | wc -l` ;
      chomp $wc_out ;
      $wc_out =~ s/\s.*$// ;
   } elsif ($fastq_fn =~ /bz2$/) {
      $wc_out = `bzcat $fastq_fn | wc -l` ;
      chomp $wc_out ;
      $wc_out =~ s/\s.*$// ;
   } else {
      $wc_out = `wc -l $fastq_fn` ;
      chomp $wc_out ;
      $wc_out =~ s/\s.*$// ;
   }

   my $num_reads = $wc_out / 4 ;
   print "$num_reads\n" ;

}
