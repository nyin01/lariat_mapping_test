
import sys
from pyfaidx import Fasta

comp_nts = {'A':'T', 'C':'G', 'T':'A', 'G':'C', 'N':'N'}
def reverse_complement(seq):
	return ''.join([comp_nts[seq[i]] for i in range(len(seq)-1,-1,-1)])

def filter_fivep_reads(fivep_to_reads, read_file, trimmed_out_path, info_out_path):

	read_sites, site_coords = {}, {}
	with open(fivep_to_reads) as fivep_file:
		for line in fivep_file:
			alignment_info = line.strip().split('\t')
			fivep_site, flag, rid, reference_start, _, read_cig = alignment_info[:6]
			for alignment_tag in alignment_info[11:]:
				if alignment_tag[:2] == 'XM':
					num_mismatch = int(alignment_tag.split(':')[-1])
			if num_mismatch == 0 and read_cig == '20M':
				reference_start = int(reference_start)-1
				bit_flags = bin(int(flag))
				is_reverse = True if len(bit_flags)>=7 and bit_flags[-5]==1 else False
				fivep_site = fivep_site[:-3]
				if rid not in read_sites:
					read_sites[rid] = set()
					site_coords[rid] = {}
				read_sites[rid].add(fivep_site)
				site_coords[rid][fivep_site] = (reference_start, reference_start+20, is_reverse)

	read_fasta = Fasta(read_file, as_raw=True)
	with open(trimmed_out_path, 'w') as trimmed_out_file, open(info_out_path, 'w') as info_out_file:
		info_out_file.write('read_id\tread_seq\tfivep_seq\tfivep_sites\tfivep_first\tread_fivep_start\tread_fivep_end\n')
		for rid in read_sites:
			fivep_info = [site_coords[rid][fp] for fp in site_coords[rid]]
			if len(set(fivep_info)) == 1:
				fivep_start, fivep_end, is_reverse = fivep_info[0]
				read_seq = read_fasta[rid][:]
				if not is_reverse:
					trim_seq, fivep_seq = read_seq[:fivep_start], read_seq[fivep_start:fivep_end]
				else:
					trim_seq, fivep_seq = read_seq[fivep_end:], reverse_complement(read_seq[fivep_start:fivep_end])
				if len(trim_seq) >= 20:
					trimmed_out_file.write('>{}\n{}\n'.format(rid, trim_seq))
					fivep_sites = ','.join(site_coords[rid].keys())
					info_out_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(rid, read_seq, fivep_seq, fivep_sites, is_reverse, fivep_start, fivep_end))
					

if __name__ == '__main__' :
	
	read_file, fivep_to_reads, trimmed_out_path, info_out_path = sys.argv[1:]
	
	filter_fivep_reads(fivep_to_reads, read_file, trimmed_out_path, info_out_path)



			

