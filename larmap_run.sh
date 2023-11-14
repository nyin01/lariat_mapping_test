#!/bin/bash

usage() {
    echo ""
    echo "Options:"
    echo "  -d, --fastq_dir <fastq_dir>               Directory containing the FASTQ files"
    echo "  -1, --read_1_file <read_1_file>           Read one (R1) FASTQ file"
    echo "  -2, --read_2_file <read_2_file>           Read two (R2) FASTQ file"
    echo "  -o, --output_dir <output_dir>             Directory for output files"
    echo "  -e, --output_base_name <output_base_name> Prefix to add to output files"
    echo "  -c, --num_cpus <num_cpus>                 Number of CPUs available"
    echo "  -i, --ref_b2index <ref_b2index>           Bowtie2 index of the full reference genome"
    echo "  -f, --ref_fasta <ref_fasta>               FASTA file of the full reference genome"
    echo "  -g, --ref_gtf <ref_gtf>                   GTF file with gene annotation of the reference genome"
    echo "  -5, --ref_5p_fasta <ref_5p_fasta>         FASTA file with sequences of first 20nt from reference 5' splice sites (first 20nt of introns)"
    echo "  -u, --ref_5p_upstream <ref_5p_upstream>   Custom file of sequences in 5nt window upstream of 5' splice sites"
    echo "  -3, --ref_3p_b2index <ref_3p_b2index>     Bowtie2 index file of last 250nt from reference 3' splice sites (last 250nt of introns)"
    echo "  -l, --ref_3p_lengths <ref_3p_lengths>     Custom file with the lengths of the sequences in ref_3p_b2index (some introns are <250nt)"
    echo "  -n, --ref_introns <ref_introns>           BED file of all introns in the reference genome"
    echo "  -m, --ref_repeatmasker <ref_repeatmasker> BED file of repetitive elements from repeatmasker"
    echo ""
    exit 1
}

exit_abnormal() {                         
    usage
    exit 1
}


# https://stackoverflow.com/questions/402377/using-getopts-to-process-long-and-short-command-line-options
echo ""
while getopts :d:1:2:o:e:c:i:f:g:5:u:3:l:n:m:-: opt; do        
    case $opt in                    
        d) 
            fastq_dir=$OPTARG
            echo "fastq_dir: $fastq_dir" ;;
        1) 
            read_one_file=$OPTARG 
            echo "read_1_file: $read_one_file" ;;
        2) 
            read_two_file=$OPTARG
            echo "read_two_file: $read_two_file" ;;
        o) 
            output_dir=$OPTARG 
            echo "output_dir: $output_dir" ;;
        e) 
            output_base_name=$OPTARG 
            echo "output_base_name: $output_base_name" ;;
        c) 
            num_cpus=$OPTARG 
            echo "num_cpus: $num_cpus" ;;
        i) 
            ref_b2index=$OPTARG 
            echo "ref_b2index: $ref_b2index" ;;
        f) 
            ref_fasta=$OPTARG
            echo "ref_fasta: $ref_fasta" ;; 
        g) 
            ref_gtf=$OPTARG 
            echo "ref_gtf: $ref_gtf" ;;
        5) 
            ref_5p_fasta=$OPTARG 
            echo "ref_5p_fasta: $ref_5p_fasta" ;;
        u) 
            ref_5p_upstream=$OPTARG 
            echo "ref_5p_upstream: $ref_5p_upstream" ;;
        3) 
            ref_3p_b2index=$OPTARG 
            echo "ref_3p_b2index: $ref_3p_b2index" ;;
        l) 
            ref_3p_lengths=$OPTARG 
            echo "ref_3p_lengths: $ref_3p_lengths" ;;
        n) 
            ref_introns=$OPTARG 
            echo "ref_introns: $ref_introns" ;;
        m) 
            ref_repeatmasker=$OPTARG 
            echo "ref_repeatmasker: $ref_repeatmasker" ;;
        -) 
            case "${OPTARG}" in
                fastq_dir)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    fastq_dir=$val
                    echo "fastq_dir: $fastq_dir" ;;
                read_1_file)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    read_one_file=$val  
                    echo "read_1_file: $read_one_file" ;;
                read_2_file)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    read_two_file=$val  
                    echo "read_2_file: $read_two_file" ;;
                output_dir)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_dir=$val  
                    echo "output_dir: $output_dir" ;;
                output_base_name)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_base_name=$val  
                    echo "output_base_name: $output_base_name" ;;
                num_cpus)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    num_cpus=$val  
                    echo "num_cpus: $num_cpus" ;;
                ref_b2index)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_b2index=$val  
                    echo "ref_b2index: $ref_b2index" ;;
                ref_fasta)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_fasta=$val 
                    echo "ref_fasta: $ref_fasta" ;;
                ref_gtf)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_gtf=$val 
                    echo "ref_gtf: $ref_gtf" ;;
                ref_5p_fasta)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_5p_fasta=$val 
                    echo "ref_5p_fasta: $ref_5p_fasta" ;;
                ref_5p_upstream)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_5p_upstream=$val 
                    echo "ref_5p_upstream: $ref_5p_upstream" ;;
                ref_3p_b2index)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_3p_b2index=$val 
                    echo "ref_3p_b2index: $ref_3p_b2index" ;;
                ref_3p_lengths)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_3p_lengths=$val 
                    echo "ref_3p_lengths: $ref_3p_lengths" ;;
                ref_introns)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_introns=$val 
                    echo "ref_introns: $ref_introns" ;;
                ref_repeatmasker)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    ref_repeatmasker=$val 
                    echo "ref_repeatmasker: $ref_repeatmasker" ;;
                *)
                    echo "Error: unrecognized option --${OPTARG}."
                    exit_abnormal   ;;
            esac;;
        :)  
            echo ""                                 
            echo "Error: -${OPTARG} requires an argument."
            exit_abnormal   ;;
        *)  
            echo ""
            echo "Error: unrecognized option -${OPTARG}."                              
            exit_abnormal   ;;
    esac
done


# Check if all required arguments are provided
if [[ -z $fastq_dir || -z $read_one_file || -z $read_two_file || -z $output_dir || -z $output_base_name || -z $num_cpus || -z $ref_b2index || -z $ref_fasta || -z $ref_gtf || -z $ref_5p_fasta || -z $ref_5p_upstream || -z $ref_3p_b2index || -z $ref_3p_lengths || -z $ref_introns || -z $ref_repeatmasker ]]; then
  echo ""
  echo "All arguments are required."
  exit_abnormal
fi

echo ""

printf "$(date +'%m/%d/%y - %H:%M:%S') | Starting lariat mapping run...\n"
script_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $script_dir
#=============================================================================#
#                                    MAPPING                                  #
#=============================================================================#

echo ""
# prepares the directories and scripts for the lariat mapping run
printf "$(date +'%m/%d/%y - %H:%M:%S') | Preparing directories...\n"
output_dir=$output_dir/$output_base_name"_lariat_mapping"
mkdir -p $output_dir/$output_base_name"_R1"
mkdir -p $output_dir/$output_base_name"_R2"

echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Processing read one file...\n"
SECONDS=0
scripts/map_lariats.sh $fastq_dir/$read_one_file \
    $output_dir/$output_base_name"_R1/"$output_base_name"_R1" \
    $num_cpus $ref_b2index $ref_fasta $ref_gtf \
    $ref_5p_fasta $ref_5p_upstream \
    $ref_3p_b2index $ref_3p_lengths
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo ""
    printf "Error: Failed to execute map_lariats.sh on read one file. Exit code: $exit_code"
    exit $exit_code
fi

echo ""
printf "$(date +'%m/%d/%y - %H:%M:%S') | Processing read two file...\n"
SECONDS=0
scripts/map_lariats.sh $fastq_dir/$read_two_file \
    $output_dir/$output_base_name"_R2/"$output_base_name"_R2" \
    $num_cpus $ref_b2index $ref_fasta $ref_gtf \
    $ref_5p_fasta $ref_5p_upstream \
    $ref_3p_b2index $ref_3p_lengths
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo ""
    printf "Error: Failed to execute map_lariats.sh on read two file. Exit code: $exit_code"
    exit $exit_code
fi

printf "$(date +'%m/%d/%y - %H:%M:%S') | Filtering results...\n"
# combines the mapping results from each sample's read one and read two files and performs post-mapping filtering before outputting the final lariat mapping results
python -u scripts/filter_lariats.py $fastq_dir $read_one_file $read_two_file \
    $output_dir $output_base_name $num_cpus $ref_b2index $ref_fasta $ref_gtf \
    $ref_5p_fasta $ref_3p_b2index $ref_3p_lengths $ref_introns $ref_repeatmasker
exit_code=$?
# Check the exit code and handle errors
if [ $exit_code -ne 0 ]; then
    echo "" 
    printf "Error: Failed to execute filter_lariats.py. Exit code: $exit_code"
    exit $exit_code
fi

printf "$(date +'%m/%d/%y - %H:%M:%S') | Finished.\n"
