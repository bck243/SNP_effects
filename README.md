# SNP_effects
Determine the effect of SNPs on amino acid sequences (e.g. are whether mutations are synonymous)

* SNPs were called on *ab initio* open reading frames (ORFs) using [bcftools](http://samtools.github.io/bcftools/bcftools.html) (mpileup, call, query) yielding a tab-delimited output with the columns: CHROM (chromosome aka ORF ID), POS (position of SNP on chromosome), REF (nucleotide identity of reference ORF at SNP location), ALT (nucleotide identity of SNP), AD (allelic depth), QUAL (quality), DP (number of reads covering the base at the call position), MQ (root mean square mapping quality). An example SNP output file is *mb1_allele_positions.tab*, included here. This tab delimited output is the [allele file] used as input for **check_ORF_frames.pl** and **get_snp_effects.pl**.

* Note that to compare SNPs across conditions downstream, read coverage should be downsampled to the levels of the least-covered condition before calling SNPs

This repository contains the following scripts: 

* **check_ORF_frames.pl**
  * perl script that checks whether a given set of ORFs are in-frame (e.g. do not start mid-codon) to determine whether they can be translated simply using *get_snp_effects.pl*. 
  * Script will break with error "Translated codon does not match reference amino acid" if ORFs are not in frame. 
  * Input: 
   * **Allele file**: a tab delimited and contains columns for ORF ID, SNP position, and reference nucleotide, in that order. Additional columns are acceptable.
   * **nucleotide fasta file**: fasta headers must be identical to ORF ID (following ">") 
   * **amino acid fasta file**: fasta headers must be identical to ORF ID (following ">") 
   * The desired output file name. 
  * Produces an output spreadsheet with columns for: ORF ID, SNP position in ORF, reference SNP nucleotide, SNP position in codon (1st, 2nd, or 3rd), codon, translated amino acid from codon, and reference amino acid from amino acid fasta file. 
  * note: will break if given a codon containing an unspecified nucleotide ("N") with error: "The reference codon is not a valid amino acid"
* **get_snp_effects.pl**
  * perl script that gives codon and amino acid translation of SNPs and reference positions, and determines whether SNPs are synonymous
  * Input: 
   * **Allele file**: a tab delimited and contains columns for ORF ID, SNP position, and reference nucleotide, in that order. Additional columns are acceptable.
   * **nucleotide fasta file**: fasta headers must be identical to ORF ID (following ">") 
   * The desired output file name. 
   * codons containing unspecified nucleotides ("N") will be translated to amino acid "X", but SNP information will still be produced
* **get_snp_effects_n.pl**
  * variant of get_snp_effects.pl that also reports nitrogeon atom flux as a result of a SNP (the difference in nitrogen atoms between the reference amino acid and new amino acid) 
