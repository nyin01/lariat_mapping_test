
from sys import argv
from os.path import join


# =============================================================================#
#                                  Constants                                  #
# =============================================================================#
HEADERS = {'fivep_info_table': "read_id	read_seq\tfivep_seq\tfivep_sites\tfivep_first\tread_fivep_start\tread_fivep_end",
           'final_info_table': "read_id\tread_seq\tchrom\tstrand\tfivep_site\tread_is_reverse\tfivep_read_start\tfivep_read_end\tthreep_site\tbp_site\tread_bp_nt\tgenomic_bp_nt\tgenomic_bp_window",
           }


# =============================================================================#
#                                  Functions                                  #
# =============================================================================#
def merge_tables(results_dir: str, results_name: str, nsplits: int) -> None:

    for table_name in HEADERS:
        with open(join(results_dir, f'{results_name}_{table_name}.txt'), 'w') as out_file:
            # write the header
            out_file.write(HEADERS[table_name] + '\n')
            # loop through the split directories
            for split_num in range(1, nsplits+1):
                # copy it into the full file
                split_name = f'{results_name}_split{split_num}'
                split_table = join(results_dir, split_name,
                                   f'{split_name}_{table_name}.txt')
                with open(split_table, 'r') as split_file:
                    next(split_file)
                    for line in split_file:
                        out_file.write(line)


def sum_readcounts(results_dir: str, results_name: str, nsplits: int) -> None:
    '''
    Get the total number of reads with 1 alignment to the reference genome for each split, add them together,
    and write the result to a .txt file in the split's results directory
    '''
    total = 0

    for split_num in range(1, nsplits+1):
        split_name = f'{results_name}_split{split_num}'
        split_readcount_file = join(
            results_dir, split_name, f'{split_name}_total_reads.txt')
        with open(split_readcount_file) as in_file:
            read_count = in_file.readline().strip().split('=')[1]
            total += int(read_count)

    if total == 0:
        raise RuntimeError(f'Failed to get total reads for {results_name}')

  # now write the total to file
    with open(join(results_dir, f'{results_name}_total_reads.txt'), 'w') as w:
        w.write(f'total_reads={total}')


# =============================================================================#
#                                    Main                                     #
# =============================================================================#
def main():

    results_dir, results_name, nsplits = argv[1:]
    nsplits = int(nsplits)

    merge_tables(results_dir, results_name, nsplits)

    sum_readcounts(results_dir, results_name, nsplits)


if __name__ == '__main__':
    main()
