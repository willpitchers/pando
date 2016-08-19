#!/usr/bin/env python


'''
Read in alignment.  Find sites that have more than one type of nucleotide
(considering only A, C, T or G).  Store their positions in variant_sites.
Add up the alignments in variant sites and write to file in fasta format.

Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm
YYYMMDD_HHMM: 20160819_1624
'''


#import modules
import argparse
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Alphabet import generic_dna



PARSER = argparse.ArgumentParser(description='Will read in an alignment in'+\
                                 ' any format, and spit out only the' +\
                                 ' variant sites to a new fasta-formatted' +\
                                 ' alignment.')
PARSER.add_argument('-f', '--filename', help='Name of input file.',
                    required=True)
PARSER.add_argument('-i', '--informat', help='File format of input.' +\
                    ' e.g., genbank, fasta, phylip', required=True)
ARGS = PARSER.parse_args()


def read_collapse(infile, informat):
    '''
    Pull out the variable sites from an alignment
    '''
    with open(infile, 'r') as input_handle:
        alignment = AlignIO.read(input_handle, informat, alphabet=generic_dna)
        summary_align = AlignInfo.SummaryInfo(alignment)
        first_seq = (alignment[0].seq)
        my_pssm = summary_align.pos_specific_score_matrix(first_seq)
        dna_bases = ['A', 'C', 'T', 'G']
        pos = 0
        variant_sites = []
        for base_dict in my_pssm:
            base_count = [base_dict[j] for j in dna_bases if base_dict[j] > 0]
            if len(base_count) > 1:
                variant_sites.append(pos)
            pos += 1
        print 'Collapsing alignment to '+str(len(variant_sites))+' sites...'
        cmd = ['alignment[:, '+str(i)+':'+str(i+1)+']' for i in variant_sites]
        snp_alignment = eval('+'.join(cmd))
        with open('core_gene_alignment_collapsed.fasta', 'w') as outfile:
            AlignIO.write(snp_alignment, outfile, 'fasta')


def main():
    '''
    Run read_collapse.
    '''
    read_collapse(ARGS.filename, ARGS.informat.replace('.', ''))


if __name__ == '__main__':
    main()
    print '\nDone.\n'
