#!/usr/local/bin/perl
#fpd150518_1948
# purpose: convert ENSEMBL chromosome names to UCSC compatible names

use strict;
use warnings;

main() ;

sub main {

   my $usage = __FILE__." [column number of chromosome name; default 1] < old_file" ;
   if ($#ARGV >= 0 && $ARGV[0] =~ /[hH]/) {
      die $usage ;
   }

   my $chrom_field = 0;
   if ($#ARGV >= 0) {$chrom_field = $ARGV[0] - 1;}

   my $ensembl2ucsc = {
      '2L'      => 'chr2L',
      '2R'      => 'chr2R',
      '3L'      => 'chr3L',
      '3R'      => 'chr3R',
      '4'       => 'chr4',
      'X'       => 'chrX',
      'Y'       => 'chrY',
      'dmel_mitochondrion_genome' => 'chrM'
   } ;

   while (my $line = <STDIN>) {
      chomp $line;
      my @t = split(/\t/, $line) ;
      if (!exists $ensembl2ucsc->{$t[$chrom_field]}) {next;}

      $t[$chrom_field] = $ensembl2ucsc->{$t[$chrom_field]} ;
      print join("\t", @t)."\n";
   }

}
