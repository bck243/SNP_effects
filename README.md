# SNP_effects
Determine the effect of SNPs on amino acid sequences (e.g. are whether mutations are synonymous)

* SNPs were called on ab initio ORFs using [bcftools](http://samtools.github.io/bcftools/bcftools.html) (mpileup, call, query) yielding a tab-delimited output with the columns: CHROM (chromosome aka ORF ID), POS (position of SNP on chromosome), REF (nucleotide identity of reference ORF at SNP location), ALT (nucleotide identity of SNP), AD (allelic depth), QUAL (quality), DP (number of reads covering the base at the call position), MQ (root mean square mapping quality). An example SNP output file is *mb1_allele_positions.tab*, included here.

* Note that to compare SNPs across conditions, read coverage should be downsampled to the levels of the least-covered condition before calling SNPs

This repository contains the following scripts: 

* **check_ORF_frames.pl**
  * perl script that checks whether a given set of ORFs are in-frame (e.g. do not start mid-codon) to determine whether they can be translated simply using *get_snp_effects.pl*. 
* **get_snp_effects.pl**
  * perl script that gives codon and amino acid translation of SNPs and reference positions, and determines whether SNPs are synonymous
