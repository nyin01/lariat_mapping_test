# Lariat-Mapping Pipeline

## Environment

This pipeline depends on the command line tools and python packages contained in the conda environment at `larmap_env` outlined in `environment.yaml`. To build the virtual environment, run `conda env create -f environment.yaml`. 

For M1 mac users: please install packages `bowtie2`, `bedtools`, and `samtools` using the command `arch -arm64 brew install [package]` before running `conda`, if any of the above pacakges has not previously been installed.

## Running the Pipeline

To run the larmap pipeline, use `./larmap_top.sh` followed by arguments listed below:

        -d, --fastq_dir                           Path to the directory containing the fastq files
        -1, --read_1_file                         Read 1 file
        -2, --read_2_file                         Read 2 file
        -e, --output_base_name                    Output base name
        -c, --num_cpus                            Number of CPUs available for the mapping run
        -i, --ref_b2index                         Path to file containing Bowtie2 index of the whole genome
        -f, --ref_fasta                           Path to Fasta file of the whole genome
        -g, --ref_gtf                             Path to GTF file providing the gene annotation of the whole genome
        -5, --ref_5p_fasta                        Path to Fasta file containing the 5' sequences of all the annotated introns (first 20bp of the introns)
        -3, --ref_3p_b2index                      Path to file containing Bowtie2 index of a custom 3' intronic genome (last 250bp of all annotated introns)
        -l, --ref_3p_lengths                      Path to file containing lengths of all the sequences in the custom 3' genome (not all introns are >=250bp, so some lengths will be shorter)
        -n, --ref_introns                         Path to BED file containing the intron annotation of the whole genome
        -m, --ref_repeatmasker                    Path to BED file containing the RepeatMasker annotation of the whole genome
        -s, --results_path                        Path to the .txt file where the final lariat mapping results will be written

For example: 
`./larmap_top.sh 
-d demo_files/demo_fastq_files_250k_bp 
-1 250k_cWT_1.fq.gz 
-2 250k_cWT_2.fq.gz 
-e WT 
-c 4 
-i demo_files/genomes/indices/bowtie2/mm39.fa 
-f demo_files/genomes/fasta_files/mm39.fa 
-g demo_files/genomes/annotations/mm39.gencode.basic.M32.annotation.gtf.gz 
-5 demo_files/reference_files/mouse/mm39.gencode.basic.M32.fivep_sites.fa 
-3 demo_files/reference_files/mouse/mm39.gencode.basic.M32.threep_sites.fa 
-l demo_files/reference_files/mouse/mm39.gencode.basic.M32.threep_seq_lens.txt 
-n demo_files/genomes/annotations/mm39.gencode.basic.M32.introns.bed.gz 
-m demo_files/genomes/annotations/mm39.repeat_masker.bed.gz 
-s demo_mapped_reads_250k.txt`

Or: `./larmap_top.sh --fastq_dir demo_files/demo_fastq_files_250k_bp ...`

A default directory `larmap_out` will be created upon running the pipeline; after completion of the pipeline it will contain three directories `larmap_out/logs`, `larmap_out/scripts`, and `larmap_out/output`, as well as a `.txt` file (as specified in the arguments) which stores the final mapping results.

## Workflow

1. The `map_lariats_top.py` script prepares the directories and scripts for the lariat mapping run. This script will read the settings and read file information from the command line arguments, and generate scripts in `larmap_out/scripts` for mapping each of the sample's read one and read two files. Each read's script is titled `larmap*[output_base_name]\_[R1/R2].sh`. In addition, it creates a `bash_all.sh` script that runs each of the bash scripts upon execution. Logs will be written in `larmap_out/logs/top_log.out`.

2. The `bash_all.sh` script runs all the individual bash scripts, which execute the mapping algorithm.

   Each `larmap*/sh` scirpt will:
   
        - Activate the `larmap_env` conda virtual environment
        - Make a directory under `larmap_out/output` named `[output_base_name]*[R1/R2]`
        - Run the `map_lariats.sh` script on the read file. This will produce three files in the read output directory: `[output_base_name]*[R1/R2]_total_reads.txt` (one line file containing count of linearly-aligned reads from read file), `[output_base_name]_[R1/R2]_fivep_info_table.txt` (intermediate file containing info on the mapping of the 5'SS sequences to the unmapped reads), and `[output_base_name]_[R1/R2]\_final_info_table.txt` (results file containing candidate lariat reads obtained after mapping the 5'SS trimmed reads to the 3'SS region sequences).
        - Output logs to its own log file within the `larmap_out/logs/child_logs` directory.

    The mapping script `map_lariats.sh` will:
   
        - Align reads to the reference genome with bowtie2 and create a BAM file of the unmapped reads with samtools
        - Convert the unmapped reads bam file to FASTA format with samtools
        - Build a bowtie2 index of the unmapped reads
        - Align a fasta file of 5'SS to the unmapped reads index
        - Trim reads with 5'SS alignments and write trimmed reads to fasta file
        - Align the trimmed reads to a Bowtie2 index of 3'SS regions
        - Take the mapped reads from step 6 and create an output file containing candidate lariat reads

3. The `merge_filter_lariats.py` script loads intron and gene information from provided annotation files, combines the mapping results from each sample's read one and read two files, and performs post-mapping filtering before outputting the final lariat mapping results. 

    The candidate lariat reads from each sample are filtered based on the following criteria: - BP is within 2bp of a splice site (likely from an intron circle, not a lariat) - 5'SS and 3'SS are not in the correct order - Read maps to a Ubiquitin gene (likely false positive due to repetitive nature of gene) - There is a valid aligment for the 3' segment upstream of the 5' segment - Both the 5'SS and the BP overlap with repetitive regions from RepeatMasker (likely false positive) - NEEDS TO BE IMPLEMENTED: Filter out template switching reads (eg. 5'SS 6bp sequence matches 6bp sequence downstream of BP) - NEEDS TO BE IMPLEMENTED: Correct BP position to account for RT skipping

    After filtering, the final set of lariat reads for all the samples are merged (so that lariats discovered in both mates are only reported once) and output to the `results_path` file. The initial columns of this tab-delimited file are any additional custom columns which were added to the info file. After these columns, the file contains the following data for each lariat read: - Gene name - Gene Ensembl ID - Gene type (e.g. protein-coding, lncRNA) - Read ID - Read sequence - Chromosome - Strand - 5' splice site coordinate - 3' splice site coordinate - Branchpoint coordinate - Read branchpoint nucleotide - Genomic branchpoint nucleotide - Genomic branchpoint sequence context (10bp window with the BP at position 5) - Branchpoint position relative to nearest 3' splice site - Total count of linearly-aligned reads for the sample the lariat read is from (useful for normalization when comparing lariat levels between samples)
