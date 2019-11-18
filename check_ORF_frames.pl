#!/usr/bin/env perl
use strict;
use warnings;
##############################
#GOALS
##############################
#fetch fasta sequences that match list of identifiers
#get nts at SNP positions for all fasta sequences
#check translate codons to aa seqs
#check that the resulting amino acids match translations done by ORF caller

##############################
#get input variables
##############################
my $allelefile = shift; #allele file (ORF ids to select; one per line)
#tab delimited; columns are orf id, SNP position, reference nucleotide, mutated nucleotide...
#(remaining columns aren't needed)
my $infa = shift; #fasta file to search through (nt sequences) e.g. "assembly.orf.ffn"
my $faa = shift; #amino acid fasta file for checking translation e.g. "aassembly.orf.faa"
my $outfile = shift;
#tab delimited output

##############################
#define remaining variables
##############################
#scalars
my $read; #1 when ORF id matches input ids
my $orf_id; #ORF id
my $snp_position; #location of snp
my $ref_snp_nt; #reference nucleotide at snp location
my $orf_id_i; #orf ID for given iteration of WHILE loop
my $snp_codon_pos; #location of SNP in codon (1st, 2nd, 3rd)
my $position1; #position of 1st character to print (of codon)
my $position2; #position of 2nd character to print (of codon)
my $position3; #position of 3rd character to print (of codon)
my $nt1;#first character to print (of codon)
my $nt2;#2nd character to print (of codon)
my $nt3;#3rd character to print (of codon)
my $codon;
my $ref_aa; #reference amino acid at snp location
my $remainder; #remainder when dividing SNP position by 3; gives codon location of SNP
my $snp_nt;
my $row_number;
my $codon_position; #position in codon (1st, 2nd, or 3rd)
#arrays
my @orf_list; #list of ORF IDs used to match fasta files
#hashes
my %orfid_fnn; #orf id, nucleotide fasta sequence
my %orfid_faa; #orf id, amino acid fasta sequence

##############################
#check input files
##############################
if($allelefile && $infa && $faa) {
} else {
  die "Usage: $0 [allele file] [nucleotide fasta file] [amino acid fasta file] [output file name]\nAllele file is tab delimited and contains columns for ORF ID, SNP position, and reference nucleotide, in that order. \n";
}
print "Reading input files \n";

##############################
#Open allele file and get ORF id hashes for
#SNP position, reference nucleotide, mutated nucleotide
##############################
open(IN, $allelefile) or die "Unable to open file $allelefile\nAllele file is tab delimited and contains columns for ORF ID, SNP position, and reference nucleotide, in that order. \n";
while(<IN>) {
    chomp;
    my ($orf_id, $snp_position, $ref_snp_nt) = split(/\t/);
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
open(IN, $infa) or die "Unable to open file $infa\n ";
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
#open aa fasta file and get orf ID/fasta hash
##############################
#loop through fasta input
open(IN, $faa) or die "Unable to open file $faa\n";
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
    $orfid_faa{$orf_id_i} = $_;
  } #end parsing fasta sequence
}
close(IN); #finish looping through nt fasta

##############################
#Open allele file again to find codons, print aa translation, and write out
##############################
open(OUT, ">".$outfile) or die "Unable to write to file $outfile\n";
#print header
print OUT join("\t", ("orf_id", "nt_snp_position", "reference_nt_snp", "snp_position_in_codon", "codon", "translated_aa", "orf_aa"))."\n";

  open(IN, $allelefile) or die "Unable to open file $allelefile\n";
  while(<IN>) {
      chomp;
      #print STDERR "before opening allele file, \$snp_position is set to".$snp_position."\n";
      my ($orf_id, $snp_position, $ref_snp_nt) = split(/\t/);
      $row_number ++;
      print STDERR "Parsing SNP ".$row_number." of ".(scalar @orf_list)."\n";
      #print STDERR "After opening allele file, \$snp_position is set to".$snp_position."\n";
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
      my $codon;
      my $ref_aa;
      if($remainder == 0){
        $snp_codon_pos = 3;
        $position1 = $snp_position -2;
        $position2 = $snp_position -1;
        $position3 = $snp_position;
        $nt1 = substr($orfid_fnn{$orf_id}, $position1 , 1 );
        $nt2 = substr($orfid_fnn{$orf_id}, $position2 , 1 );
        $nt3 = substr($orfid_fnn{$orf_id}, $position3 , 1 );
      }elsif($remainder == 1){ #remainder = 1/3
        $snp_codon_pos = 1;
        $position1 = $snp_position;
        $position2 = $snp_position +1;
        $position3 = $snp_position +2;
        $nt1 = substr($orfid_fnn{$orf_id}, $position1 , 1 );
        $nt2 = substr($orfid_fnn{$orf_id}, $position2 , 1 );
        $nt3 = substr($orfid_fnn{$orf_id}, $position3 , 1 );
      }elsif($remainder == 2){ #remainder = 2/3
        $snp_codon_pos = 2;
        $position1 = $snp_position -1;
        $position2 = $snp_position;
        $position3 = $snp_position +1;
        $nt1 = substr($orfid_fnn{$orf_id}, $position1 , 1 );
        $nt2 = substr($orfid_fnn{$orf_id}, $position2 , 1 );
        $nt3 = substr($orfid_fnn{$orf_id}, $position3 , 1 );
      }else{
        die "The remainder is not a multiple of 1/3\n";
      }
      #subtract one because indexing in perl starts at 0
      #print STDERR "The SNP is at character ".($snp_position+1)."\n";
      #add 1 back for human readable printing
      $snp_nt = substr($orfid_fnn{$orf_id}, $snp_position , 1 );
      $codon = ($nt1.$nt2.$nt3);
      #print STDERR "The SNP is ".$snp_nt."\n";
      #print STDERR "The remainder is ".$remainder."\n";
      #print STDERR "The codon is ".$codon."\n";
      #print STDERR "The nt seq is ".$orfid_fnn{$orf_id}."\n";

      ##############################
      #translate codon
      ##############################
      if($codon eq "GCT" or $codon eq "GCC" or $codon eq "GCA" or $codon eq "GCG"){
        $ref_aa = "A";
      }elsif($codon eq "CGT" or $codon eq "CGC" or $codon eq "CGA" or $codon eq "CGG" or $codon eq "AGA" or $codon eq "AGG"){
        $ref_aa = "R";
      }elsif($codon eq "AAT" or $codon eq "AAC"){
        $ref_aa = "N";
      }elsif($codon eq "GAT" or $codon eq "GAC"){
        $ref_aa = "D";
      }elsif($codon eq "TGT" or $codon eq "TGC"){
        $ref_aa = "C";
      }elsif($codon eq "CAA" or $codon eq "CAG"){
        $ref_aa = "Q";
      }elsif($codon eq "GAA" or $codon eq "GAG"){
        $ref_aa = "E";
      }elsif($codon eq "GGT" or $codon eq "GGC" or $codon eq "GGA" or $codon eq "GGG"){
        $ref_aa = "G";
      }elsif($codon eq "CAT" or $codon eq "CAC"){
        $ref_aa = "H";
      }elsif($codon eq "ATT" or $codon eq "ATC" or $codon eq "ATA"){
        $ref_aa = "I";
      }elsif($codon eq "TTA" or $codon eq "TTG" or $codon eq "CTT" or $codon eq "CTC" or $codon eq "CTA" or $codon eq "CTG"){
        $ref_aa = "L";
      }elsif($codon eq "AAA" or $codon eq "AAG"){
        $ref_aa = "K";
      }elsif($codon eq "ATG"){
        $ref_aa = "M";
      }elsif($codon eq "TTT" or $codon eq "TTC"){
        $ref_aa = "F";
      }elsif($codon eq "CCT" or $codon eq "CCC" or $codon eq "CCA" or $codon eq "CCG"){
        $ref_aa = "P";
      }elsif($codon eq "TCT" or $codon eq "TCC" or $codon eq "TCA" or $codon eq "TCG" or $codon eq "AGT" or $codon eq "AGC"){
        $ref_aa = "S";
      }elsif($codon eq "ACT" or $codon eq "ACC" or $codon eq "ACA" or $codon eq "ACG"){
        $ref_aa = "T";
      }elsif($codon eq "TGG"){
        $ref_aa = "W";
      }elsif($codon eq "TAT" or $codon eq "TAC"){
        $ref_aa = "Y";
      }elsif($codon eq "GTT" or $codon eq "GTC" or $codon eq "GTA" or $codon eq "GTG"){
        $ref_aa = "V";
      }elsif($codon eq "TAA" or $codon eq "TGA" or $codon eq "TAG"){
        $ref_aa = "STOP";
      }else{
        die "The reference codon is not a valid amino acid\n";
      }
      #print STDERR "The reference codon translates to ".$ref_aa."\n";

      ##############################
      #get aa at SNP position from aa fasta file
      ##############################
      my $snp_aa_pos;
      $snp_aa_pos = int($snp_position/3);
      #subtract one because indexing in perl starts at 0
      #print STDERR "The SNP is at amino acid position ".($snp_aa_pos+1)."\n";
      #add 1 back for human readable printing
      my $snp_aa;
      $snp_aa = substr($orfid_faa{$orf_id}, $snp_aa_pos, 1 );
      #print STDERR "The reference amino acid at the SNP location is ".$snp_aa."\n";
      #print STDERR "The amino acid seq is ".$orfid_faa{$orf_id}."\n";

      ##############################
      #break if codon translation does not match amino acid from ORF
      ##############################

      if($ref_aa ne $snp_aa){
        die "Translated codon does not match reference amino acid\n";
      }
      ##############################
      #write out
      ##############################
      #output is: "orf_id", "nt_snp_position", "nt_snp", "snp_position_in_codon", "codon", "translated_aa", "orf_aa"
      print OUT join("\t", ($orf_id, ($snp_position+1), $ref_snp_nt, $codon_position, $codon, $ref_aa, $snp_aa))."\n";
      print STDERR "Parsed SNP ".$row_number." of ".(scalar @orf_list)."\n----------------------------\n";
  }
  close(IN);
  print STDERR "Parsed allele file \n";
close(OUT);
print STDERR "Done. \n";
