#!/bin/bash

#=============================================================================#
#                                  Arguments                                  #
#=============================================================================#

# RNA-seq fastq read file
READ_FILE=$1
# Output directory with sample name prefix
OUTPUT_BASE=$2
# Number of CPUs to use
CPUS=$3
# Genome Bowtie2 index base name
GENOME_BOWTIE2_INDEX=$4
# Reference genome FASTA file
GENOME_FASTA=$5
# GTF file containing gene annotatin for mapping genome
GTF_FILE=$6
# FASTA file of 5' splice sites (first 20nts of all introns)
FIVEP_FASTA=$7
# Custom file of sequences in 5nt window upstream of 5'ss
FIVEP_UPSTREAM=$8 
# Bowtie2 index of 3' splice sites genome (last 250nts of all introns)
THREEP_BOWTIE2_INDEX=$9
# 3' splice site lengths
THREEP_LENGTHS="${10}"

#=============================================================================#
#                                    Calls                                    #
#=============================================================================#
### Map filtered reads to genome and keep unmapped reads. Lariat reads crossing the brachpoint will not be able to map to the gene they're from
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping reads and extracting unmapped reads...\n"
output_bam=$OUTPUT_BASE"_mapped_reads.bam"
unmapped_bam=$OUTPUT_BASE"_unmapped_reads.bam"
bowtie2 --end-to-end --sensitive --score-min L,0,-0.24 -k 1 --n-ceil L,0,0.05 --threads $CPUS -x $GENOME_BOWTIE2_INDEX -U $READ_FILE \
	| samtools view --bam --with-header > $output_bam
samtools view --bam --with-header --require-flags 4 $output_bam > $unmapped_bam
mapped_read_count=$(samtools view --count --exclude-flags 4 $output_bam)
echo "total_reads=$mapped_read_count" > $OUTPUT_BASE"_total_reads.txt"

### Produce files for Shapeshifter analysis
# echo ""
# echo "SHAPESHIFTER"
# # do bedtools covreage on output_bam to produce files required for shapeshifter
# sorted_bam=$OUTPUT_BASE"_mapped_reads_sorted.bam"
# echo "samtools sort"
# samtools sort $output_bam -o $sorted_bam -@ $CPUS
# coverage=$OUTPUT_BASE"_coverage.bedgraph"
# echo "bedtools coverage"
# bedtools genomecov -trackline -bg -ibam $sorted_bam > $coverage
# # Extract the first 50 nt of each intron from the input intron bed file
# echo "extract the first 50 nt"
# awk '{print $1"\t"$2"\t"$2+50"\t"$4}' $coverage > $OUTPUT_BASE"_50nt.bed"
# # Sort and remove duplicates based on the first three columns (chromosome, start, end)
# echo "sort and remove duplicates"
# sort -k1,1 -k2,2n -k3,3n -u $OUTPUT_BASE"_50nt.bed" > $OUTPUT_BASE"_50nt_sorted.bed"
# # Use bedtools intersect to get coverage for the 5' intron segments
# echo "bedtools intersect"
# bedtools intersect -wo -a $coverage -b $OUTPUT_BASE"_50nt_sorted.bed" > $OUTPUT_BASE"_intron_cov.txt"


### Create fasta file of unmapped reads 
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Creating fasta file of unmapped reads...\n"
unmapped_fasta=$OUTPUT_BASE"_unmapped_reads.fa"
samtools fasta $unmapped_bam > $unmapped_fasta
samtools faidx $unmapped_fasta

### Build a bowtie index of the unmapped reads
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Building bowtie index of unmapped fasta...\n"
bowtie2-build --large-index --threads $CPUS $unmapped_fasta $unmapped_fasta > /dev/null

### Align unmapped reads to fasta file of all 5' splice sites (first 20nts of introns)
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' splice sites to reads...\n"
fivep_to_reads=$OUTPUT_BASE"_fivep_to_reads.sam"
bowtie2 --end-to-end --sensitive --no-unal -f --k 10000 --threads $CPUS -x $unmapped_fasta -U $FIVEP_FASTA \
	| samtools view > $fivep_to_reads

### Extract reads with a mapped 5' splice site and trim it off
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Finding 5' read alignments and trimming reads...\n"
fivep_info_table=$OUTPUT_BASE"_fivep_info_table.txt"
fivep_trimmed_reads=$OUTPUT_BASE"_fivep_mapped_reads_trimmed.fa"
python scripts/filter_fivep_alignments.py $unmapped_fasta $fivep_to_reads $FIVEP_UPSTREAM $fivep_trimmed_reads $fivep_info_table

### Map 5' trimmed reads to 3' sites (last 250nts of introns)
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Mapping 5' trimmed reads to 3' sites...\n"
trimmed_reads_to_threep=$OUTPUT_BASE"_fivep_reads_trimmed_mapped_to_threep.sam"
bowtie2 --end-to-end --sensitive -k 10 --no-unal --threads $CPUS -f -x $THREEP_BOWTIE2_INDEX -U $fivep_trimmed_reads \
	| samtools view > $trimmed_reads_to_threep

### Filter 3' splice site alignments and output info table, including the branchpoint site
echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Analyzing 3' alignments and outputting lariat table...\n"
python scripts/filter_threep_alignments.py $trimmed_reads_to_threep $THREEP_LENGTHS $fivep_info_table $GTF_FILE $GENOME_FASTA $OUTPUT_BASE

### Delete all intermediate/uneeded files that were created throughout this process
wait
rm $output_bam
rm $unmapped_bam
rm $unmapped_fasta* $fivep_to_reads* $fivep_trimmed_reads $trimmed_reads_to_threep*  
