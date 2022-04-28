# Populus-simonii-mitogenome

This repository contains some Perl scripts for phylogenetic analysis of mitochondrial genomes.

1. Step1_Get_gene_matrix.pl 

   Usage: perl Step1_Get_gene_matrix.pl -i Species_name.txt -o Gene_matrix.txt
   
   The Species_name.txt file should contain all the file name of your selected data (Fasta format).
   
   [Species_name.txt](https://github.com/bichangwei/Populus-simonii-mitogenome/files/8579230/Species_name.txt)

   The purpose of this script is to generate the following gene matrix.  
   
   [Gene_matrix.txt](https://github.com/bichangwei/Populus-simonii-mitogenome/files/8579182/Gene_matrix.txt)
   
   Then, the users need to select the conserved single genes from the gene matrix, and output the conserved gene list for the next step.
   
   [ConservedGene.txt](https://github.com/bichangwei/Populus-simonii-mitogenome/files/8579206/ConservedGene.txt)

2. Step2_Generate_the_phylogeny_fa.pl

   Usage: perl Step2_Generate_the_phylogeny_fa.pl -g ConservedGene.txt -o Prefix
   
   This script will generate two main results:
   
   (1) One file named $prefix.phylogeny.fa in Fasta format for constructing the phylogenetic tree.
    
   (2) IQTREE results. The parameters for IQTREE are: iqtree -s $output.muscle.alignment.afa -m MFP -B 1000 --bnni -T AUTO
   
3. All the data can be found in the Example directory.
   
   
   
