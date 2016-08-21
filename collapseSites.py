#!/usr/bin/env python


'''
Read in alignment.  Find sites that have more than one type of nucleotide
(considering only A, C, T or G).  Store their positions in variant_sites.
Add up the alignments in variant sites and write to file in fasta format.

Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm
YYYMMDD_HHMM: 20160821_1923
'''


#import modules
import argparse
from Bio import AlignIO
from Bio import SeqIO
import sys
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from multiprocessing import Pool




PARSER = argparse.ArgumentParser(description='Will read in an alignment in'+\
                                 ' any format, and spit out only the' +\
                                 ' variant sites to a new fasta-formatted' +\
                                 ' alignment.')
PARSER.add_argument('-f', '--filename', help='Name of input file.',
                    required=True)
PARSER.add_argument('-i', '--informat', help='File format of input.' +\
                    ' e.g., genbank, fasta, phylip', required=True)
PARSER.add_argument('-t', '--threads', help='Number of threads.' +\
                    ' default=72', default=72, type=int, required=True)

ARGS = PARSER.parse_args()

def summary_aln(aln_params):
    '''
    Take an alignment and coordinates (unpacked from aln_params).  Take a 
    sub-alignment bounded by the coordinates.  Calculate a pssm for the 
    sub-alignment.  Find the variable sites using the pssm.  Take a 
    sub-sub-alignment for each variable site.  Bind all the sub-sub-alignments
    to return a sub-alignment of just variable sites (SNVs).
    '''
    (aln, i, j) = aln_params
    alignment = aln[:, i:j]
    summary_align = AlignInfo.SummaryInfo(alignment)
    first_seq = (alignment[0].seq)
    my_pssm = summary_align.pos_specific_score_matrix(first_seq)
    dna_bases = ['A', 'C', 'T', 'G']
    pos = 0
    variant_sites = []
    for base_dict in my_pssm:
        base_count = [base_dict[k] for k in dna_bases if base_dict[k] > 0]
        if len(base_count) > 1:
            variant_sites.append(pos)
        pos += 1
    variant_cols = []
    for s in variant_sites:
        variant_cols.append(alignment[:,s:s+1])
    aln = MultipleSeqAlignment(variant_cols[0])
    for t in variant_cols[1:]:
        aln += t
    return aln


def read_collapse(infile, informat):
    '''
    Pull out the variable sites from an alignment.  Do this in parallel. 
    '''
    with open(infile, 'r') as input_handle:
        alignment = AlignIO.read(input_handle, informat, alphabet=generic_dna)
        print 'Read '+infile+' for determining variant positions...'
        print 'Converted sequences to uppercase...'
        for record in alignment:
            record.seq = record.seq.upper()
        aln_len = alignment.get_alignment_length()
        window = aln_len//(ARGS.threads-1)
        aln_ranges = []
        c = 0
        while c+window < aln_len:
            aln_ranges.append((alignment, c, c+window))
            c += window+1
        else:
            aln_range = [c, aln_len]
            aln_ranges.append((alignment, c, aln_len))
        print 'Processing the alignment in parallel to find variant sites.'
        p = Pool(ARGS.threads)
        results = p.map(summary_aln, aln_ranges)
        aln = MultipleSeqAlignment(results[0])
        for i in results:
            aln += i
        with open('core_gene_alignment_collapsed.fasta', 'w') as outfile:
            AlignIO.write(aln, outfile, 'fasta')
            print 'Written collapsed alignment to' +\
                  ' core_gene_alignment_collapsed.fasta'


def main():
    '''
    Run read_collapse.
    '''
    read_collapse(ARGS.filename, ARGS.informat.replace('.', ''))


if __name__ == '__main__':
    main()
    print '\nDone.\n'
