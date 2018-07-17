#!/usr/local/bin/perl
#fpd 150518_2144
# Purpose: extract transcript info from BDGP (+ERCC + INTACT) GTF File

use strict;
use warnings;
main() ;

sub main {

   print join("\t", qw/transcript_id gene_id transcript_name gene_name gene_biotype/)."\n";
   while (my $line = <STDIN>) {
      chomp $line;
      if ($line !~ /\ttranscript\t/ &&
          $line !~ /\tERCC/ &&
          $line !~ /\tINTACT/) {next;}
      my @t = split(/\t/, $line) ;
      my ($gene_id) =  ($line =~ /gene_id \"([^"]+)\"/) ;
      my ($tx_id) =  ($line =~ /transcript_id \"([^"]+)\"/) ;
      my ($gene_name) = ($line =~ /gene_name \"([^"]+)\"/) ;
      my ($tx_name) = ($line =~ /transcript_name \"([^"]+)\"/) ;

      my ($gene_biotype) = ($line =~ /gene_biotype \"([^"]+)\"/) ;

      if (!defined $gene_name) {$gene_name = $gene_id ;}
      if (!defined $tx_name) {$tx_name = $tx_id ;}
      if (!defined $gene_biotype) {$gene_biotype = 'unk'; }

      if ($gene_id =~ /^ERCC/) {$tx_id = $gene_id;} 
      print join("\t", $tx_id, $gene_id, $tx_name, $gene_name, $gene_biotype)."\n";
   }

}

#3R	FlyBase	transcript	722370	722621	.	-	.	gene_id "FBgn0085804"; gene_version "1"; transcript_id "FBtr0114258"; transcript_version "1"; gene_name "CR41571"; gene_source "FlyBase"; gene_biotype "pseudogene"; transcript_name "CR41571-RA"; transcript_source "FlyBase"; transcript_biotype "pseudogene";
#3R	FlyBase	transcript	835381	2503907	.	+	.	gene_id "FBgn0267431"; gene_version "1"; transcript_id "FBtr0346770"; transcript_version "1"; gene_name "CG45784"; gene_source "FlyBase"; gene_biotype "protein_coding"; transcript_name "CG45784-RA"; transcript_source "FlyBase"; transcript_biotype "protein_coding";
