=head1 NAME

extract_interpro_from_flybase_GFF.pl

=head1 VERSION

fpd160607_1244

=head1 AUTHOR

Fred P. Davis (fredpdavis@gmail.com)

=head1 SYNOPSIS

extract_interpro_from_flybase_GFF.pl makes a table of InterPro annotations

=cut


#2L	FlyBase	gene	1143201	1147395	.	+	.	ID=FBgn0031307;Name=MFS3;fullname=Major Facilitator Superfamily Transporter 3;Alias=CG4726;Ontology_term=SO:0000010,SO:0000087,GO:0005316,GO:0016021,GO:0055085;Dbxref=FlyBase:FBan0004726,FlyBase_Annotation_IDs:CG4726,GB_protein:AAF51413,GB:AY095062,GB_protein:AAM11390,GB:BI162517,GB:BI162518,UniProt/TrEMBL:Q9VPX2,INTERPRO:IPR011701,INTERPRO:IPR016196,OrthoDB7_Drosophila:EOG7MDH8M,OrthoDB7_Diptera:EOG7SZ6KQ,EntrezGene:33292,INTERPRO:IPR020846,BDGP_clone:RE01809,OrthoDB7_Insecta:EOG7SJR43,OrthoDB7_Arthropoda:EOG7HXR52,OrthoDB7_Metazoa:EOG789C9Z,FlyAtlas:CG4726-RA,Fly-FISH:CG4726,GenomeRNAi:33292;gbunit=AE014134;derived_computed_cyto=21F1-21F1

use strict;
use warnings;

main() ;

sub main {

   if ($#ARGV >= 0) {die "USAGE: ".__FILE__." < flybase.GFF" ;}

   print "gene_id\tdomain\n";
   while (my $line = <STDIN>) {
      chomp $line;
      if ($line !~ /\tgene\t/ || $line !~ /INTERPRO/) {next;}
      my ($fbgn) = ($line =~ /\tID=(FBgn[0-9]+);/) ;
      my (@iprdoms) = ($line =~ /INTERPRO:(IPR[0-9]+)/g) ;
      foreach my $i ( 0 .. $#iprdoms) {
         print $fbgn."\t".$iprdoms[$i]."\n"; }
   }

}
