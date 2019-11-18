#!/usr/bin/env perl
use strict;
use warnings;
##############################
#GOALS
##############################
#previously checked that ORFs (assembly.orf.ffn) are in-frame and can be translated directly
#here, translated ORFs to get reference amino acid at SNP position
#use SNP identity to get mutatated amino acid
#determine whether SNPs are synonymous

##############################
#define variables for input files
##############################
my $allelefile = shift; #allele file (ORF ids to select; one per line)
#tab delimited; columns are orf id, SNP position, reference nucleotide, mutated nucleotide...
#(remaining columns aren't needed)
my $infa = shift; #fasta file to search through (nt sequences)  e.g. "assembly.orf.ffn"
my $outfile = shift; # e.g."translated_snps.tab";
#tab delimited output

##############################
#define remaining variables
##############################
#scalars
my $read; #1 when ORF id matches input ids
my $orf_id; #ORF id
my $snp_position; #location of snp
my $ref_snp_nt; #reference nucleotide at snp location
my $new_snp_nt; #new nucleotide at snp location
my $orf_id_i; #orf ID for given iteration of WHILE loop
my $ref_aa; #reference amino acid at snp location
my $new_aa; #new amino acid at snp location
my $remainder; #remainder when dividing SNP position by 3; gives codon location of SNP
my $nt1;#first character to print (of reference codon)
my $nt2;#2nd character to print (of reference codon)
my $nt3;#3rd character to print (of reference codon)
my $snp_codon_pos; #location of SNP in codon (1st, 2nd, 3rd)
my $ref_codon; #reference codon
my $new_codon; #reference codon after mutation
my $position1; #position of 1st character to print (of reference codon)
my $position2; #position of 2nd character to print (of reference codon)
my $position3; #position of 3rd character to print (of reference codon)
my $codon_position; #position of SNP in codon (1st, 2nd, or 3rd)
my $is_synonymous; #0 = non-synonymous; 1 = synonymous
my $row_number;
#arrays
my @orf_list;
#hashes
my %orfid_fnn; #orf id, nucleotide fasta sequence
my %orfid_faa; #orf id, amino acid fasta sequence

##############################
#check input files
##############################
if($allelefile && $infa) {
} else {
  die "Usage: $0 [allele file] [nucleotide fasta file] [output file name]\nAllele file is tab delimited and contains columns for ORF ID, SNP position, and reference nucleotide, in that order. \n";
}
print "Reading input files \n";

##############################
#Open allele file and get ORF id hashes for
#SNP position, reference nucleotide, mutated nucleotide
##############################
open(IN, $allelefile) or die "Unable to open file $allelefile\nAllele file is tab delimited and contains columns for ORF ID, SNP position, and reference nucleotide, in that order. \n";
while(<IN>) {
    chomp;
    my ($orf_id, $snp_position, $ref_snp_nt, $new_snp_nt) = split(/\t/);
    #my $pattern = $orf_id.$snp_position;
    #pattern must be ORF + postion because there are redundant ORF IDs
    push(@orf_list, $orf_id);
}
close(IN);
print STDERR "Obtained ORF ID's \n";
print STDERR (scalar @orf_list)." SNPs detected \n----------------------------\n";

##############################
#open nt fasta file and get orf ID/fasta hash
##############################
#loop through fasta input
open(IN, $infa) or die "Unable to open file $infa\n";
#loop through rows of fasta file
while(<IN>) {
  #look for >; check if it matches orf id
  if(/^>/) { #if the line starts with a >
    $read = 0; #set "read" to zero (eg skip printing)
    foreach my $orf_id (@orf_list) { #check if ID matches
        if(/$orf_id/) { #if it does;
            $read = 1; #set read to 1 (line will be printed while $read =1)
            $orf_id_i = $orf_id;
        }
    }
  }
  #if fasta is a match, print it
  if($read) {
    #print STDERR $_; #orf id; eg ">contig_501235_1_1467_+"
  }
  #get sequence line output; substring to proper nt position
  if($read & !/^>/) { #if $read = 1 (ID matched) AND there's no carrot (fasta only)
    $orfid_fnn{$orf_id_i} = $_;
    #print STDERR $orfid_fnn{$orf_id_i};
  } #end parsing fasta sequence
}
close(IN); #finish looping through nt fasta

##############################
#Open allele file again to find codons, print aa translation, and write out
##############################

open(OUT, ">".$outfile) or die "Unable to write to file $outfile\n";
#print header
print OUT join("\t", ("orf_id", "snp_position", "reference_nucleotide", "first_new_nucleotide", "reference_codon", "new_codon", "snp_position_in_codon","translated_reference_codon", "translated_new_codon", "synonymous"))."\n";

  open(IN, $allelefile) or die "Unable to open file $allelefile\n";
  print STDERR "Translating SNP-containing codons \n";
  while(<IN>) {
      chomp;
      my ($orf_id, $snp_position, $ref_snp_nt, $new_snp_nt) = split(/\t/);
      $row_number ++;
      print STDERR "Parsing ORF ".$row_number." of ".(scalar @orf_list)."\n";
      ##############################
      #Get SNP position
      ##############################
      $snp_position = $snp_position -1;
      my $remainder;
      $remainder = ($snp_position+1) % 3;
      my $codon_position; #position in codon (1st, 2nd, or 3rd)
      if ($remainder == 0){
        $codon_position = 3;
      }elsif ($remainder ==1){
        $codon_position = 1;
      }elsif ($remainder == 2){
        $codon_position = 2;
      }else{
        die "Impossible codon position\n";
      }
      my $snp_codon_pos; #location of SNP in codon (1st, 2nd, 3rd)
      my $position1; #position of 1st character to print (of codon)
      my $position2; #position of 2nd character to print (of codon)
      my $position3; #position of 3rd character to print (of codon)
      my $nt1;#first character to print (of codon)
      my $nt2;#2nd character to print (of codon)
      my $nt3;#3rd character to print (of codon)
      my $ref_codon;
      my $new_codon;
      my $ref_aa;
      my $new_aa;
      if($remainder == 0){
        $snp_codon_pos = 3;
        $position1 = $snp_position -2;
        $position2 = $snp_position -1;
        $position3 = $snp_position;
        $nt1 = substr($orfid_fnn{$orf_id}, $position1 , 1 );
        $nt2 = substr($orfid_fnn{$orf_id}, $position2 , 1 );
        $nt3 = substr($orfid_fnn{$orf_id}, $position3 , 1 );
        $new_codon = ($nt1.$nt2.$new_snp_nt);
      }elsif($remainder == 1){ #remainder = 1/3
        $snp_codon_pos = 1;
        $position1 = $snp_position;
        $position2 = $snp_position +1;
        $position3 = $snp_position +2;
        $nt1 = substr($orfid_fnn{$orf_id}, $position1 , 1 );
        $nt2 = substr($orfid_fnn{$orf_id}, $position2 , 1 );
        $nt3 = substr($orfid_fnn{$orf_id}, $position3 , 1 );
        $new_codon = ($new_snp_nt.$nt2.$nt3);
      }elsif($remainder == 2){ #remainder = 2/3
        $snp_codon_pos = 2;
        $position1 = $snp_position -1;
        $position2 = $snp_position;
        $position3 = $snp_position +1;
        $nt1 = substr($orfid_fnn{$orf_id}, $position1 , 1 );
        $nt2 = substr($orfid_fnn{$orf_id}, $position2 , 1 );
        $nt3 = substr($orfid_fnn{$orf_id}, $position3 , 1 );
        $new_codon = ($nt1.$new_snp_nt.$nt3);
      }else{
        die "The remainder is not a multiple of 1/3\n";
      }
      #subtract one because indexing in perl starts at 0
      #print STDERR "The SNP is at character ".($snp_position+1)."\n";
      #add 1 back for human readable printing
      $ref_snp_nt = substr($orfid_fnn{$orf_id}, $snp_position , 1 );
      $ref_codon = ($nt1.$nt2.$nt3);
      #print STDERR "The SNP is ".$ref_snp_nt."\n";
      #print STDERR "The remainder is ".$remainder."\n";
      #print STDERR "The codon is ".$ref_codon."\n";
      #print STDERR "The nt seq is ".$orfid_fnn{$orf_id}."\n";

      ##############################
      #translate reference codon
      ##############################
      if($ref_codon eq "GCT" or $ref_codon eq "GCC" or $ref_codon eq "GCA" or $ref_codon eq "GCG"){
        $ref_aa = "A";
      }elsif($ref_codon eq "CGT" or $ref_codon eq "CGC" or $ref_codon eq "CGA" or $ref_codon eq "CGG" or $ref_codon eq "AGA" or $ref_codon eq "AGG"){
        $ref_aa = "R";
      }elsif($ref_codon eq "AAT" or $ref_codon eq "AAC"){
        $ref_aa = "N";
      }elsif($ref_codon eq "GAT" or $ref_codon eq "GAC"){
        $ref_aa = "D";
      }elsif($ref_codon eq "TGT" or $ref_codon eq "TGC"){
        $ref_aa = "C";
      }elsif($ref_codon eq "CAA" or $ref_codon eq "CAG"){
        $ref_aa = "Q";
      }elsif($ref_codon eq "GAA" or $ref_codon eq "GAG"){
        $ref_aa = "E";
      }elsif($ref_codon eq "GGT" or $ref_codon eq "GGC" or $ref_codon eq "GGA" or $ref_codon eq "GGG"){
        $ref_aa = "G";
      }elsif($ref_codon eq "CAT" or $ref_codon eq "CAC"){
        $ref_aa = "H";
      }elsif($ref_codon eq "ATT" or $ref_codon eq "ATC" or $ref_codon eq "ATA"){
        $ref_aa = "I";
      }elsif($ref_codon eq "TTA" or $ref_codon eq "TTG" or $ref_codon eq "CTT" or $ref_codon eq "CTC" or $ref_codon eq "CTA" or $ref_codon eq "CTG"){
        $ref_aa = "L";
      }elsif($ref_codon eq "AAA" or $ref_codon eq "AAG"){
        $ref_aa = "K";
      }elsif($ref_codon eq "ATG"){
        $ref_aa = "M";
      }elsif($ref_codon eq "TTT" or $ref_codon eq "TTC"){
        $ref_aa = "F";
      }elsif($ref_codon eq "CCT" or $ref_codon eq "CCC" or $ref_codon eq "CCA" or $ref_codon eq "CCG"){
        $ref_aa = "P";
      }elsif($ref_codon eq "TCT" or $ref_codon eq "TCC" or $ref_codon eq "TCA" or $ref_codon eq "TCG" or $ref_codon eq "AGT" or $ref_codon eq "AGC"){
        $ref_aa = "S";
      }elsif($ref_codon eq "ACT" or $ref_codon eq "ACC" or $ref_codon eq "ACA" or $ref_codon eq "ACG"){
        $ref_aa = "T";
      }elsif($ref_codon eq "TGG"){
        $ref_aa = "W";
      }elsif($ref_codon eq "TAT" or $ref_codon eq "TAC"){
        $ref_aa = "Y";
      }elsif($ref_codon eq "GTT" or $ref_codon eq "GTC" or $ref_codon eq "GTA" or $ref_codon eq "GTG"){
        $ref_aa = "V";
      }elsif($ref_codon eq "TAA" or $ref_codon eq "TGA" or $ref_codon eq "TAG"){
        $ref_aa = "STOP";
      }else{
        die "The reference codon is not a valid amino acid\n";
      }
      print STDERR "The reference codon translates to ".$ref_aa."\n";

      ##############################
      #translate new codon
      ##############################
      if($new_codon eq "GCT" or $new_codon eq "GCC" or $new_codon eq "GCA" or $new_codon eq "GCG"){
        $new_aa = "A";
      }elsif($new_codon eq "CGT" or $new_codon eq "CGC" or $new_codon eq "CGA" or $new_codon eq "CGG" or $new_codon eq "AGA" or $new_codon eq "AGG"){
        $new_aa = "R";
      }elsif($new_codon eq "AAT" or $new_codon eq "AAC"){
        $new_aa = "N";
      }elsif($new_codon eq "GAT" or $new_codon eq "GAC"){
        $new_aa = "D";
      }elsif($new_codon eq "TGT" or $new_codon eq "TGC"){
        $new_aa = "C";
      }elsif($new_codon eq "CAA" or $new_codon eq "CAG"){
        $new_aa = "Q";
      }elsif($new_codon eq "GAA" or $new_codon eq "GAG"){
        $new_aa = "E";
      }elsif($new_codon eq "GGT" or $new_codon eq "GGC" or $new_codon eq "GGA" or $new_codon eq "GGG"){
        $new_aa = "G";
      }elsif($new_codon eq "CAT" or $new_codon eq "CAC"){
        $new_aa = "H";
      }elsif($new_codon eq "ATT" or $new_codon eq "ATC" or $new_codon eq "ATA"){
        $new_aa = "I";
      }elsif($new_codon eq "TTA" or $new_codon eq "TTG" or $new_codon eq "CTT" or $new_codon eq "CTC" or $new_codon eq "CTA" or $new_codon eq "CTG"){
        $new_aa = "L";
      }elsif($new_codon eq "AAA" or $new_codon eq "AAG"){
        $new_aa = "K";
      }elsif($new_codon eq "ATG"){
        $new_aa = "M";
      }elsif($new_codon eq "TTT" or $new_codon eq "TTC"){
        $new_aa = "F";
      }elsif($new_codon eq "CCT" or $new_codon eq "CCC" or $new_codon eq "CCA" or $new_codon eq "CCG"){
        $new_aa = "P";
      }elsif($new_codon eq "TCT" or $new_codon eq "TCC" or $new_codon eq "TCA" or $new_codon eq "TCG" or $new_codon eq "AGT" or $new_codon eq "AGC"){
        $new_aa = "S";
      }elsif($new_codon eq "ACT" or $new_codon eq "ACC" or $new_codon eq "ACA" or $new_codon eq "ACG"){
        $new_aa = "T";
      }elsif($new_codon eq "TGG"){
        $new_aa = "W";
      }elsif($new_codon eq "TAT" or $new_codon eq "TAC"){
        $new_aa = "Y";
      }elsif($new_codon eq "GTT" or $new_codon eq "GTC" or $new_codon eq "GTA" or $new_codon eq "GTG"){
        $new_aa = "V";
      }elsif($new_codon eq "TAA" or $new_codon eq "TGA" or $new_codon eq "TAG"){
        $new_aa = "STOP";
      }else{
        print STDERR "The new codon is not a valid amino acid\n";
      }
      print "The new codon translates to ".$new_aa."\n";

      ##############################
      #check if new amino acid is synonymous
      ##############################
      my $is_synonymous;
      if ($new_aa eq $ref_aa){
        $is_synonymous = "yes";
      }else{
        $is_synonymous = "no";
      }

      ##############################
      #write out
      ##############################
      #output is: "orf_id", "nt_snp_position", "nt_snp", "remainder", "codon", "translated_aa", "orf_aa"
      #"orf_id", "snp_position", "reference_nucleotide", "first_new_nucleotide", "reference_codon", "new_codon", "snp_position_in_codon","translated_reference_codon", "translated_new_codon", "synonymous"
      print OUT join("\t", ($orf_id, ($snp_position+1), $ref_snp_nt, $new_snp_nt, $ref_codon, $new_codon, $codon_position, $ref_aa, $new_aa, $is_synonymous))."\n";
      print STDERR "Parsed ORF ".$row_number." of ".(scalar @orf_list)."\n----------------------------\n";
  }
  close(IN);
  print STDERR "Parsed allele file \n";
close(OUT);
print STDERR "Done. \n";
