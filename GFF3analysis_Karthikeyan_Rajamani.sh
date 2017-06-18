#!/bin/bash

#Description -comparing genome architecture between C.elegans, D.melanogaster and H.sapiens .gff3 files; output is produced in Results.txt
#Author Karthikeyan Rajamani
#Usage ./GFF3analysis_Karthikeyan_Rajamani.sh 
#download .gff3 files

curl -o C.elegans.gff3.gz ftp://ftp.ensemblgenomes.org/pub/metazoa/release-34/gff3/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.34.gff3.gz
curl -o D.melanogaster.gff3.gz ftp://ftp.ensemblgenomes.org/pub/metazoa/release-34/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.34.chr.gff3.gz
curl -o H.sapiens.gff3.gz ftp://ftp.ensembl.org/pub/release-87/gff3/homo_sapiens/Homo_sapiens.GRCh38.87.chr.gff3.gz
gunzip *.gff3.* 
# Header for the files
echo "Processing......"
awk 'BEGIN{print "Species\tTotalGenes\tTotalTranscripts\tTranscriptsPerGene\tGenesWith1Transcript\tGenesWith2OrMoreTranscripts\tMaxNumberOfTranscriptsPerGene\tGenomeSize(Mb)\tGeneDensity(genes/Mb)"}' >Results.txt
for file in *.gff3; do
#get species name
species_name=$(basename $file .gff3)
#get total genes
awk '$1!="dmel_mitochondrion_genome"&&$1!="MtDNA"&&$1!="MT" {print $0}' $file>mitocondrial_eliminate.txt
egrep 'protein_coding' mitocondrial_eliminate.txt >protien_coding_only.txt
gene_count=$(awk '$3=="gene" {sum=sum+1}END{print sum}' protien_coding_only.txt)
#get total transcripts
transcript_count=$(awk '$3=="mRNA" {sum=sum+1}END{print sum}' protien_coding_only.txt)
#get transcripts/gene
transcripts_per_gene=$(awk -v gc=$gene_count -v tc=$transcript_count 'BEGIN {print tc/gc}')
#5. Total number of genes with 1 transcript
awk '$3=="mRNA" {print $9}' protien_coding_only.txt>nuclear_genes_9thcol.txt
awk 'BEGIN{OFS="\t"}{print gensub(/.*transcript_id=([^;]*).*/,"\\1",1,$0),gensub(/.*Parent=gene:([^;]*).*/,"\\1",1,$0)}' nuclear_genes_9thcol.txt >GeneID_ProteinID.txt
cut -f2 GeneID_ProteinID.txt | uniq -c >Gene_Transcript_count.txt
single_transcript=$(awk '$1=="1"{sum=sum+1}END{print sum}' Gene_Transcript_count.txt)
#6. Total number of genes with >1 transcripts
greaterthanone_transcript=$(awk '$1>1{sum=sum+1}END{print sum}' Gene_Transcript_count.txt)
#7. Maximum number of transcripts per gene
Max_Transcripts=$(sort -n -k1 Gene_Transcript_count.txt|tail -n1|awk '{print $1}')
#8. Size of NUCLEAR genome in megabases (1Mb=1000000 bp)
Nuc_Size_MB=$(awk '$1=="##sequence-region" && $2!="dmel_mitochondrion_genome" {sum=sum+ $4} END{print sum/1000000}' $file)
#9. Gene density (genes/Mb)
Genedensity=$(awk -v gc=$gene_count -v mb=$Nuc_Size_MB  'BEGIN{print gc/mb}')
awk -v v1=$species_name -v v2=$gene_count -v v3=$transcript_count -v v4=$transcripts_per_gene -v v5=$single_transcript -v v6=$greaterthanone_transcript -v v7=$Max_Transcripts -v v8=$Nuc_Size_MB -v v9=$Genedensity 'BEGIN {OFS="\t";print v1,v2,v3,v4,v5,v6,v7,v8,v9}' >>Results.txt
rm mitocondrial_eliminate.txt protien_coding_only.txt nuclear_genes_9thcol.txt GeneID_ProteinID.txt Gene_Transcript_count.txt 
done
echo "Done-View your Results.txt"
	
	
	
	