#!/usr/local/bin/perl
#fpd 150721_0833
# Purpose: make intron and exon BED files from BDGP (+ERCC + INTACT) GTF File

use strict;
use warnings;
use File::Temp qw/tempfile/ ;
main() ;

sub main {

   my $usage = __FILE__." gtf_fn output_prefix" ;

   my $specs = {} ;
   $specs->{gtf_fn}     = $ARGV[0] || die $usage ;
   $specs->{outprefix}  = $ARGV[1] || die $usage ;
   $specs->{bedtools_bin} = "~/bin/bedtools-v2.15.0/bedtools" ;

   $specs->{exonbed_fn} = $specs->{outprefix}.".gene_exons.bed" ;
   $specs->{intronbed_fn} = $specs->{outprefix}.".gene_introns.bed" ;


# first 3 lines of BDGP GTF:
#      3R      FlyBase gene    722370  722621  .       -       .       gene_id "FBgn0085804"; gene_version "1"; gene_name "CR41571"; gene_source "FlyBase"; gene_biotype "pseudogene";
#      3R      FlyBase transcript      722370  722621  .       -       .       gene_id "FBgn0085804"; gene_version "1"; transcript_id "FBtr0114258"; transcript_version "1"; gene_name "CR41571"; gene_source "FlyBase"; gene_biotype "pseudogene"; transcript_name "CR41571-RA"; transcript_source "FlyBase"; transcript_biotype "pseudogene";
#      3R      FlyBase exon    722370  722621  .       -       .       gene_id "FBgn0085804"; gene_version "1"; transcript_id "FBtr0114258"; transcript_version "1"; exon_number "1"; gene_name "CR41571"; gene_source "FlyBase"; gene_biotype "pseudogene"; transcript_name "CR41571-RA"; transcript_source "FlyBase"; transcript_biotype "pseudogene"; exon_id "FBtr0114258-E1"; exon_version "1";

# 1. make Genes and Exons BED file
   my $temp_fh = {} ;
   my $temp_fn = {} ;

   my $proper_chrom = {
      "2L" => 1,
      "2R" => 1,
      "3L" => 1,
      "3R" => 1,
      "4" => 1,
      "X" => 1,
      "Y" => 1,
      "dmel_mitochondrion_genome" => 1,
   } ;

   ($temp_fh->{genebed}, $temp_fn->{genebed}) = tempfile() ;
   ($temp_fh->{exonbed}, $temp_fn->{exonbed}) = tempfile() ;
   ($temp_fh->{genebed_sort}, $temp_fn->{genebed_sort}) = tempfile() ;
      close($temp_fh->{genebed_sort});
   ($temp_fh->{exonbed_sort}, $temp_fn->{exonbed_sort}) = tempfile() ;
      close($temp_fh->{exonbed_sort});
   ($temp_fh->{exonbed_sort_merged}, $temp_fn->{exonbed_sort_merged}) = tempfile() ;
      close($temp_fh->{exonbed_sort_merged}) ;
   ($temp_fh->{intronbed}, $temp_fn->{intronbed}) = tempfile() ;
      close($temp_fh->{intronbed});

   my $gene2exon ;
   open(GTF, $specs->{gtf_fn}) ;
   while (my $line = <GTF>) {

      if ($line =~ /^#/)                                { next; }

      my @t = split(/\t/, $line) ;
      if (($t[2] ne 'gene' && $t[2] ne 'exon') ||
          !exists $proper_chrom->{$t[0]})               { next; }

      my ($gene_id) =  ($line =~ /gene_id \"([^"]+)\"/) ;


      my $out_fh = $temp_fh->{genebed};
      if ($t[2] eq 'exon') {
         $out_fh = $temp_fh->{exonbed};
         push @{$gene2exon->{$gene_id}}, $t[0]."\t".($t[3] - 1)."\t".$t[4] ;
      }

      print {$out_fh} join("\t", $t[0], ($t[3] - 1), $t[4], $gene_id)."\n";

   }
   close($temp_fh->{genebed}) ;
   close($temp_fh->{exonbed}) ;

#   die $temp_fn->{genebed}." ".$temp_fn->{exonbed} ;


# 2. Sort BED files

   my $tcom = "sort -k1,1 -k2,2n -k3,3n ".$temp_fn->{exonbed}." > ".
              $temp_fn->{exonbed_sort} ;
   system($tcom) ;

   $tcom = "sort -k1,1 -k2,2n -k3,3n ".$temp_fn->{genebed}." > ".
              $temp_fn->{genebed_sort} ;
   system($tcom) ;


# 3. subtract to get intron BED file

   print STDERR "Computing intron BED file\n" ;
   $tcom = $specs->{bedtools_bin}." subtract -a ".$temp_fn->{genebed_sort}.
              " -b ".$temp_fn->{exonbed_sort}." | sort -k1,1 -k2,2n -k3,3n > ".
              $temp_fn->{intronbed} ;
   system($tcom) ;

   system("cp ".$temp_fn->{intronbed}." ".$specs->{intronbed_fn}) ;


# 4. Merge exons belonging to same gene
   my ($tfh, $tfn) = tempfile() ; close($tfh) ;
   print STDERR "Merging exons for gene #:  0" ;
   my $n = 1 ;
   foreach my $gene (sort keys %{$gene2exon}) {

      print STDERR "\b"x(length($n - 1)).$n ;

      open($tfh, ">$tfn") ;
      print {$tfh} join("\n", @{$gene2exon->{$gene}})."\n";
      close($tfh) ;

      my $tcom = "sort -k1,1 -k2,2n -k3,3n $tfn | ".
                 $specs->{bedtools_bin}." merge | ".
                 "sed 's/\$/	$gene/' >> ".$temp_fn->{exonbed_sort_merged} ;
      system($tcom) ;
      $n++ ;
   }
   print STDERR "X\n";

   $tcom = "sort -k1,1 -k2,2n -k3,3n ".$temp_fn->{exonbed_sort_merged}.
           " > ".$specs->{exonbed_fn} ;
   system($tcom) ;

# 5. cleanup!
   foreach my $x (keys %{$temp_fh}) { unlink $temp_fh->{$x} ; }

}
