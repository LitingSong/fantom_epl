#!/bin/bash

cd /home1/GENE_proc/SONGLITING/FANTOM/bed_files
 
# gencode annotation file: http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
# 1. positive strand
cat /home1/songlt/ref_genome/genome/gencode_GRCh38.p13/gencode.v38.annotation.gtf |awk -F '\t' '{if ($3=="gene" && $7=="+" ){print $0}}' |awk -F ';' '{if ($2==" gene_type \"protein_coding\""){print $0}}' > pos_protein_coding.txt 

## get TSS
cat pos_protein_coding.txt |cut -f1,4,7 > pos_protein_coding.TSS.txt

## get genesymbol
cat pos_protein_coding.txt |cut -d '"' -f6 > pos_protein_coding.symbol.txt

## combine
paste pos_protein_coding.TSS.txt pos_protein_coding.symbol.txt > pos_protein_coding.comb.txt


# 2. negative strand
cat /home1/songlt/ref_genome/genome/gencode_GRCh38.p13/gencode.v38.annotation.gtf |awk -F '\t' '{if ($3=="gene" && $7=="-" ){print $0}}' |awk -F ';' '{if ($2==" gene_type \"protein_coding\""){print $0}}' > neg_protein_coding.txt       

## get TSS
cat neg_protein_coding.txt |cut -f1,5,7 > neg_protein_coding.TSS.txt

## get genesymbol
cat neg_protein_coding.txt |cut -d '"' -f6 > neg_protein_coding.symbol.txt

## combine
paste neg_protein_coding.TSS.txt neg_protein_coding.symbol.txt > neg_protein_coding.comb.txt


# 3. promoter.bed
cat pos_protein_coding.comb.txt neg_protein_coding.comb.txt > protein_coding.comb.txt
cat protein_coding.comb.txt| awk -v OFS="\t" '{print $1,$2-500,$2+500,$4,0,$3,0,0,0,0,0,0}' > protein_coding.promoter.bed

rm ./*txt
