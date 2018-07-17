#!/usr/local/bin/perl
#fpd 160106_2334

# tally A, C, G, T and transcript length for all entries in fasta file

use strict;
use warnings;
main() ;

sub main {

    my $cur_seq = {} ;
    print join("\t", ("transcript_id", "tx_length",
                      "perc_a","perc_c","perc_g","perc_t"))."\n";
    while (my $line = <STDIN>) {
      chomp $line;
      if ($line =~ /^>/) {
         if (exists $cur_seq->{name}) {
            print_info($cur_seq) ;
            map {delete $cur_seq->{$_}; } keys %{$cur_seq} ;
         }
         my $name = $line;
         $name =~ s/^>// ;
         $name =~ s/ .*// ;
         $cur_seq->{name} = $name ;
      } else {
         $cur_seq->{seq} .= $line ;
      }
    }
    print_info($cur_seq) ;

}

sub print_info {
   my $in = shift ;

   my $seq_length = length($in->{seq}) ;

   my $nts = {};
   foreach my $i (0 .. ($seq_length - 1)) {
      $nts->{substr($in->{seq},$i,1)}++ ; }

   foreach my $nt (keys %{$nts}) {
      $nts->{$nt} /= $seq_length ;
      $nts->{$nt} = sprintf("%.2f", 100 * $nts->{$nt}) ;
   }

   print join("\t", $in->{name}, $seq_length,
              $nts->{A}, $nts->{C}, $nts->{G}, $nts->{T})."\n";
                  
}
