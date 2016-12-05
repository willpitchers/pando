#!/usr/bin/env python


'''
Run requests on a list of isolate IDs.
Read in metadata from Excel spreadsheet.
Add MALDI-ToF results, submit lab, and sp.
Return a supermatrix of metadata and a tree.
Specific for MDU folder structures and QC
(but could be adapted for non-MDU folder structures).
Optionally run roary analysis.
Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm
YYYMMDD_HHMM: 20161012_2300

Acknowledgements:
Torsten Seemann (pando relies heavily on Torsten's tools)
Anders Goncalves da Silva (advice on multiprocessing)
Dieter Bulach (advice on setting up fripan on the MDU server)
Jason Kwong (roary2fripan, advised to run Kraken on the contigs)
Susan Ballard (ongoing feature requests)
Tim Stinear (advised to use more-than-kraken to bin isolates into species)
Andrew Page (for all roary tools)
Ben Howden (feature requests)

Needs Jason Kwong's roary2fripan.py script in the working dir.

In summary, move into a directory where you want to run this script.
Enter the following two commands at the prompt:

wget raw.githubusercontent.com/MDU-PHL/pando/master/pando.py
wget raw.githubusercontent.com/MDU-PHL/pando/master/collapseSites.py
wget raw.githubusercontent.com/kwongj/roary2fripan/master/roary2fripan.py

After completing the above, (on an MDU server) now run pando using:
time nice python pando.py -i isos.txt
'''


import os
import argparse
import sys
import string
import random
from subprocess import Popen, PIPE
import shlex
import tempfile
import shutil
import itertools
import glob
from collections import defaultdict
from multiprocessing import Pool
from Bio import AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import pandas as pd
from ete3 import Tree


VERSION = 'pando version 2.3.4'


# set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Run exploratory analyses.')
PARSER.add_argument('-i', '--mdu_read_IDs', help="One MDU-ID per line\
                    (required). Put in same folder as run folder.", #any path?
                    required=True)
PARSER.add_argument('-n', '--new_IDs', help='Enter IDs (space delimited) that\
                    you wish to \'flag-if-new\' in the final table.',
                    nargs='+', required=False)
PARSER.add_argument('-w', '--wgs_qc', help='Path to WGS\
                    QC. Default=\'/mnt/seq/MDU/QC/\'',
                    default='/mnt/seq/MDU/QC/', required=False)
PARSER.add_argument('-d', '--delete_tempdirs', help='Delete tempdirs created\
                    during run? Default=\'yes\'.', default='yes',
                    required=False)
PARSER.add_argument("-t", "--threads", help='Number of threads,\
                    default=\'72\'', default=72, type=int, required=False)
PARSER.add_argument('-a', '--andi_run', help='Run andi phylogenomic analysis?\
                    Default=\'no\'', default='no', required=False)
PARSER.add_argument('-r', '--roary_run', help='Run roary pangenome analysis?\
                    Default=\'no\'', default='no', required=False)
PARSER.add_argument('-m', '--metadata_run', help='Gather metadata for all\
                    isolates? Default=\'no\'', default='no', required=False)
PARSER.add_argument('-s', '--model_andi_distance', help='Substitution model.\
                    \'Raw\', \'JC\', or \'Kimura\'. Default=\'JC\'.',
                    default='JC', required=False)
PARSER.add_argument('-c', '--percent_cutoff', help='For abricate, call the\
                    gene \'present\' if greater than this value and \'maybe\'\
                    if less than this value. Default=95. NB: 100 percent=\
                    \'100\', not \'1\'.', default=95, type=int, required=False)
PARSER.add_argument('-x', '--excel_spreadsheet', help='Parse excel spreadsheet\
                    of metadata (.xlsx format)to extract LIMS data.\
                    The data must start on line 5 (1-based indexing, as per\
                    Excel line numbers), which is default for LIMS, and\
                    contain columns with the labels \'MDU sample ID\',\
                    \'Species identification (MALDI-TOF)\',\
                    \'Species identification (Subm. lab)\' and \'Submitter\'.',
                    required=False)

ARGS = PARSER.parse_args()


if ARGS.threads > 72:
    print 'Number of requested threads must be less than 72. Exiting now.'
    sys.exit()
print '\nStarting '+VERSION+'...'
print str(ARGS.threads) +' CPU processes requested.'


#Add MLST schemes to force their usage if that species is encountered
#Only force schemes if there are two (e.g., A baumannii and E coli)
FORCE_MLST_SCHEME = {"Acinetobacter baumannii": "abaumannii",
                     #"Citrobacter freundii": "cfreundii",
                     #"Cronobacter": "cronobacter",
                     #"Enterobacter cloacae": "ecloacae",
                     "Escherichia coli": "ecoli",
                     #"Klebsiella oxytoca": "koxytoca",
                     #"Klebsiella pneumoniae": "kpneumoniae",
                     #"Pseudomonas aeruginosa": "paeruginosa"
                    }


class Isolate(object):
    '''
    An Isolate class takes the sequence ID and associates the data with that
    ID.

        attributes: ID
    '''

    def __init__(self, ID):
        '''
        From the QC path, initialise the object with an MDU-ID.
        '''
        self.ID = ID

    def assembly(self):
        '''
        Store the path to the assembly, or tell the user if it doesn't exist.
        Write missing IDs to file. Print the filename of missing isolates.
        '''
        isolate_qc_contigs = ARGS.wgs_qc+self.ID+'/contigs.fa'
        if os.path.exists(isolate_qc_contigs):
            return isolate_qc_contigs
        else:
            print isolate_qc_contigs+' does not exist'

    def assembly_metrics(self):
        '''
        Return the contig metrics by running 'fa -t'.
        '''
        os.system('fa -t '+ARGS.wgs_qc+self.ID+'/contigs.fa > '+self.ID+\
                  '_metrics.txt')
        metrics = [line.rstrip().split('\t') for line in open(self.ID+\
                   '_metrics.txt').readlines()]
        metrics = [i[1:] for i in metrics]
        metrics[0] = ['metricsContigs_'+i for i in metrics[0]]
        metrics = dict(zip(metrics[0], metrics[1]))
        os.system('rm '+self.ID+'_metrics.txt')
        metrics_df = pd.DataFrame([metrics], index=[self.ID])
        return metrics_df

    def get_yield(self):
        '''
        Return the read metrics stored in yield.tab QC file.
        '''
        yield_data = [line.rstrip().split('\t') for line in open(ARGS.wgs_qc+\
                      self.ID+'/yield.tab').readlines()]
        yield_data = dict(('metricsReads_'+i[0], i[1]) for i in yield_data[1:])
        yield_data_df = pd.DataFrame([yield_data], index=[self.ID])
        return yield_data_df

    def reads(self):
        '''
        Store the path to the QCd (trimmomatic-ed and FLASHed) reads.
        '''
        pass

    def abricate(self):
        '''
        Store the path to abricate results.
        '''
        abricate_path = ARGS.wgs_qc+self.ID+'/abricate.tab'
        if os.path.exists(abricate_path):
            ab_data = pd.read_table(abricate_path, sep='\t', header=0)
        else: #run abricate.
            os.system('abricate '+ARGS.wgs_qc+self.ID+'/contigs.fa > '+abricate_path)
            ab_data = pd.read_table(abricate_path, sep='\t', header=0)
        nrows = ab_data.shape[0]
        genes = ab_data['GENE'].tolist()
        cov = ab_data['%COVERAGE'].tolist()
        yes = []
        maybe = []
        for i in range(0, nrows):
            if cov[i] >= ARGS.percent_cutoff:
                yes.append('resgene_'+genes[i])
            else:
                maybe.append('resgene_'+genes[i])
        y = {key:'yes' for (key) in yes}
        m = {key: 'maybe' for (key) in maybe}
        #Join dictionaries m and y to form ab_results
        ab_results = dict(itertools.chain(y.iteritems(), m.iteritems()))
        #Convert to pandas dataframe
        abricate_results = pd.DataFrame([ab_results], index=[self.ID])
        return abricate_results

    def mlst(self, species, assembly):
        '''
        Store the MLST. If the species is listed in FORCE_MLST_SCHEME redo
        MLST.  If the mlst.tab exists and the species is not in FORCE, use
        the existing mlst.tab.
        '''
        if species in FORCE_MLST_SCHEME:
            cmd = 'mlst --scheme '+FORCE_MLST_SCHEME[species]+' --quiet '+\
                   assembly
            args_mlst = shlex.split(cmd)
            #should write a proc function to do this call as it's used often
            proc = Popen(args_mlst, stdout=PIPE)
            output = proc.stdout.read()
            mlst = output.rstrip().split('\n')
            mlst = [line.strip().split('\t') for line in filter(None, mlst)]
            header = mlst[0]
            data = mlst[1]
            ncol = len(header)
            mlst_formatted_dict = {'MLST_Scheme': data[1], 'MLST_ST': data[2]}
            k = 1
            for i in range(3, ncol):
                locus = header[i]+'('+data[i]+')'
                mlst_formatted_dict['MLST_Locus'+str(k)] = locus
                k += 1
        else:
            mlst_tab = ARGS.wgs_qc+self.ID+'/mlst.tab'
            if os.path.exists(mlst_tab):
                mlst = [line.rstrip().split('\t') for line in
                        open(mlst_tab).readlines()][0]
                mlst_formatted_dict = {'MLST_Scheme': mlst[1],
                                       'MLST_ST': mlst[2]}
                k = 1
                for i in range(3, len(mlst)):
                    mlst_formatted_dict['MLST_Locus'+str(k)] = mlst[i]
                    k += 1
            else:
                cmd = 'mlst --quiet '+assembly
                args_mlst = shlex.split(cmd)
                proc = Popen(args_mlst, stdout=PIPE)
                output = proc.stdout.read()
                out = output.rstrip().split('\t')[1:]
                ncol = len(out)
                mlst_formatted_dict = {'MLST_Scheme': out[0],
                                       'MLST_ST': out[1]}
                k = 1
                for i in range(3, ncol):
                    mlst_formatted_dict['MLST_Locus'+str(k)] = out[i]
                    k += 1
        mlst_results = pd.DataFrame([mlst_formatted_dict], index=[self.ID])
        return mlst_results

    def kraken_reads(self):
        '''
        Get the kraken best hit from reads.
        '''
        #Pipe these commands together
        cmd_grep = "grep -P '\tS\t' "+ARGS.wgs_qc+'/'+self.ID+"/kraken.tab"
        cmd_sort = 'sort -k 1 -r'
        cmd_head = 'head -3'

        #Split the cmds using shlex, store in args
        args_grep = shlex.split(cmd_grep)
        args_sort = shlex.split(cmd_sort)
        args_head = shlex.split(cmd_head)

        #Pipe the output of one args to another
        proc1 = Popen(args_grep, stdout=PIPE)
        proc2 = Popen(args_sort, stdin=proc1.stdout, stdout=PIPE)
        proc3 = Popen(args_head, stdin=proc2.stdout, stdout=PIPE)

        output = proc3.stdout.read()
        kraken = output.rstrip().split('\n')
        kraken = [line.strip().split('\t') for line in filter(None, kraken)]
        return kraken

    def kraken_contigs(self):
        '''
        Get the kraken best hit from assemblies.
        '''
        #Pipe these commands together
        cmd_kraken = 'nice kraken --threads 2 --db /bio/db/kraken/minikraken'+\
                     ' --fasta-input '+ARGS.wgs_qc+'/'+self.ID+'/contigs.fa'
        cmd_krk_r = 'kraken-report'
        cmd_grep = "grep -P '\tS\t'"
        cmd_sort = 'sort -k 1 -r'
        cmd_head = 'head -3'

        #Split the cmds using shlex, store in args
        args_kraken = shlex.split(cmd_kraken)
        args_krk_report = shlex.split(cmd_krk_r)
        args_grep = shlex.split(cmd_grep)
        args_sort = shlex.split(cmd_sort)
        args_head = shlex.split(cmd_head)

        #Pipe the output of one args to another
        proc1 = Popen(args_kraken, stdout=PIPE)
        proc2 = Popen(args_krk_report, stdin=proc1.stdout, stdout=PIPE,
                      stderr=PIPE)
        proc3 = Popen(args_grep, stdin=proc2.stdout, stdout=PIPE, stderr=PIPE)
        proc4 = Popen(args_sort, stdin=proc3.stdout, stdout=PIPE, stderr=PIPE)
        proc5 = Popen(args_head, stdin=proc4.stdout, stdout=PIPE, stderr=PIPE)

        output = proc5.stdout.read()
        kraken = output.rstrip().split('\n')
        kraken = [line.strip().split('\t') for line in filter(None, kraken)]
        return kraken

    def prokka_contigs(self):
        '''
        Returns the command for running prokka on the isolate.
        '''
        cmd = 'prokka --centre X --compliant --locustag '+self.ID+\
              ' --prefix '+self.ID+' --fast --quiet --outdir %s/'+self.ID+\
              ' --cpus 2 --norrna --notrna --force '+ self.assembly()
        return cmd


def shortened_ID():
    '''
    Create a random 9 character (alphanumeric) tag for use as a tempfile
    name.
    '''
    chars = string.ascii_uppercase + string.digits
    return ''.join(random.SystemRandom().choice(chars) for _ in range(10))

def make_tempdir():
    '''
    Make a temporary directory.
    '''
    tempdir = tempfile.mkdtemp(dir='.')
    return tempdir

def read_file_lines(file_name):
    '''
    Read in lines of a file, store in a list.
    '''
    file_contents = [line.rstrip().split() for line in
                     open(file_name).readlines()]
    return file_contents

def lower_tri(full_matrix):
    '''
    Take a symmetrical matrix, convert it to a lower triangle _Matrix object.
    '''
    lower_triangle = []
    names = []
    k = 2
    for i in full_matrix:
        lower_triangle.append(map(float, i[1:k]))
        names.append(i[0])
        k += 1
    matrix = _DistanceMatrix(names, lower_triangle)
    return matrix

def get_isolate_request_IDs(ID_file):
    '''
    Reads in the MDU IDs from the request IDs file and returns IDs as a list.
    ID file must contain only one ID per line.
    '''
    IDs = list(set(filter(None, [ID.rstrip() for ID in open(ID_file, 'r').readlines()])))
    return IDs

def new_IDs(IDs):
    '''
    Find the new isolate IDs to flag in the final supermatrix of metadata.'
    '''
    new_isos_to_flag = []
    for ID in IDs:
        new_isos = [iso for iso in glob.glob(ARGS.wgs_qc+ID+'*')]
        if len(new_isos) > 0:
            for j in new_isos:
                new_isos_to_flag.append(j)
    new_isos_to_flag = [i.split('/')[-1] for i in new_isos_to_flag]
    return new_isos_to_flag

def isolates_available(IDs):
    '''
    Of the request IDs, find which ones are actually available for analysis.
    '''
    avail = []
    not_avail = []
    for ID in IDs:
        isos = [iso for iso in glob.glob(ARGS.wgs_qc+ID+'*')]
        if len(isos) > 0:
            for j in isos:
                avail.append(j)
        else:
            not_avail.append(ARGS.wgs_qc+ID)
    #Use this set function to remove any duplicates
    avail = sorted(list(set(avail)))
    print '\nFound folders with the following IDs:'
    print '\n'.join(avail)
    not_avail = sorted(list(set(not_avail)))
    if len(not_avail) > 0:
        outf = os.path.splitext(ARGS.mdu_read_IDs)[0]+'_not_found.txt'
        with open(outf, 'w') as not_found:
            not_found.write('\n'.join(not_avail))
        print '\nCould not find IDs:\n', '\n'.join(not_avail)+'\n'
    return avail

def kraken_contigs_multiprocessing(iso):
    '''
    Perform a kraken search on a contigs file using .kraken_contigs().
    Using the output of that search (a 2D array), return a pandas dataframe
    of results.
    '''
    #This dict creation step could/should be done in the class?
    ID = Isolate(iso)
    krak_cntgs = ID.kraken_contigs()
    krk_df = kraken_results_df_creator(krak_cntgs, 'cntgs')
    kraken_df = pd.DataFrame([krk_df], index=[iso])
    return kraken_df

def kraken_reads_multiprocessing(iso):
    '''
    Retrieve results from a kraken search on reads (file in QC folder) using
    .kraken_reads().  Using the output of that search (a 2D array), return a
    pandas dataframe of results.
    '''
    #This dict creation step could/should be done in the class?
    ID = Isolate(iso)
    krak_reads = ID.kraken_reads()
    krk_df = kraken_results_df_creator(krak_reads, 'reads')
    kraken_df = pd.DataFrame([krk_df], index=[iso])
    return kraken_df

def kraken_results_df_creator(kraken_hits, rds_or_cntgs):
    '''
    Take the 2D array from a kraken search and return result as a dictionary.
    '''
    #This dict creation step could/should be done in the class?
    dict_hits = {}
    k = 1
    for i in range(0, len(kraken_hits)):
        dict_hits['sp_krkn'+str(k)+'_'+rds_or_cntgs] =\
                  kraken_hits[i][5].lstrip()
        dict_hits['sp_krkn'+str(k)+'_'+rds_or_cntgs+'_pc'] = kraken_hits[i][0]
        k += 1
    return dict_hits

def abricate_multiprocessing(iso):
    '''
    Return a pd.Dataframe of abricate results stored in abricate.tab.
    '''
    ID = Isolate(iso)
    return ID.abricate()

def metricsContigs_multiprocessing(iso):
    '''
    Return a pd.Dataframe of contig metrics after running 'fa -t'.
    '''
    ID = Isolate(iso)
    return ID.assembly_metrics()

def metricsReads_multiprocessing(iso):
    '''
    Return a pd.Dataframe of read metrics stored in yield.tab.
    '''
    ID = Isolate(iso)
    return ID.get_yield()

def excel_metadata(xlsx_file):
    '''
    Read in an excel spreadsheet.
    '''
    xlsx = pd.read_excel(xlsx_file, skiprows=4, index_col=0)
    print 'Excel spreadsheet:\n'
    print xlsx
    print ''
    return xlsx

def prokka(params):
    '''
    Run prokka on the isolate. Unpack params to get iso and assembly_tempdir.
    '''
    (x, y) = params
    ID = Isolate(x)
    cmd = ID.prokka_contigs() % y
    print cmd
    os.system(cmd)

def roary(base, sp, gffs):
    ''''
    Run roary on the gff files output by prokka.
    '''
    cmd = 'nice roary -f '+base+'_'+sp+'_roary -e -n -v -z -p '+\
          str(ARGS.threads)+' '+gffs
    print '\nRunning roary for '+sp+' with the command:\n'+cmd
    os.system(cmd)

def pw_calc(aln_seq_coords):
    '''
    Calculate the absolute distance between two sequences, in a range of seqs.
    '''
    df_dist = []
    for m in aln_seq_coords:
        (aln, i, j) = m
        print i, j
        name = aln[i].id
        x = len([y for y in zip(aln[i].seq,
                aln[j].seq) if len(set(y)) > 1 and 'N' not in set(y) and '-' not in set(y) and '?' not in set(y)])
        df = pd.DataFrame([{aln[j].id: x}], index=[name])
        df_dist.append(df)
    return pd.concat(df_dist, axis=1)


def main():
    '''
    Read in the MDU-IDs from file. For each ID, instantiate an object of
    class Isolate.  This class associates QC data with the ID tag.
    Move the contigs for all isolates into a tempdir, with a temp 9-character
    filename.  Run andi phylogenomics on all the contig sets.  Infer an NJ tree
    using Bio Phylo from the andi-calculated distance matrix.  Correct the
    negative branch lengths in the NJ tree using ETE3.  Export the tree to
    file. Gather and combine the metadata for each ID as a super-matrix.
    Optionally, add LIMS metadata to the super-matrix from a LIMS excel
    spreadsheet option (adds MALDI-ToF, Submitting Lab ID, Submitting Lab
    species guess) and/or use the flag-if-new to highlight
    'new' isolates.  Export the tree and metadata to .csv, .tsv/.tab file.
    Export the 'isolates not found' to text file too.
    '''

    #i) read in the IDs from file
    IDs = get_isolate_request_IDs(ARGS.mdu_read_IDs)
    #base should be a global, given that it is used in other functions too.
    base = os.path.splitext(ARGS.mdu_read_IDs)[0]

    #ii) Return a folder path to the QC data for each available ID
    #    using a wildcard search of the ID in IDs in ARGS.wgs_qc path.
    iso_paths = isolates_available(IDs)
    #Drop the path and keep the folder name
    isos = [i.split('/')[-1] for i in iso_paths]

    #iii) make tempdir to store the temp_contigs there for 'andi' analysis.
    assembly_tempdir = make_tempdir()

    #iv) flag new isolate IDs.
    new_ids = None
    if ARGS.new_IDs != None:
        new_ids = new_IDs(ARGS.new_IDs)
    #v) read in the LIMS metadata
    if ARGS.excel_spreadsheet != None:
        xls_table = excel_metadata(ARGS.excel_spreadsheet)

    #vi) Copy contigs to become temp_contigs into tempdir, only if andi
    #requested.
    #Translation dict to store {random 9-character filename: original filename}
    iso_ID_trans = {}
    #Dict to store each isolate under each consensus species#####maybe delete
    isos_grouped_by_cons_spp = defaultdict(list)
    for iso in isos:
        #Instantiate an Isolate class for each isolate in isolates
        sample = Isolate(iso)
        #Next, we could just use iso_path+/contigs.fa, but that would skip
        #the if os.path.exists() test in sample.assembly(iso).
        assembly_path = sample.assembly()
        short_id = shortened_ID()
        #Store key,value as original_name,short_id for later retrieval.
        iso_ID_trans[iso] = short_id
        if 'y' in ARGS.andi_run.lower():
            cmd = 'ln -s '+assembly_path+' '+assembly_tempdir+'/'+short_id+\
                  '_contigs.fa'
            os.system(cmd)
            print 'Creating symlink:', cmd
    with open(base+'_temp_names.txt', 'w') as tmp_names:
        print '\nTranslated isolate IDs:\nShort\tOriginal'
        for key, value in iso_ID_trans.items():
            print value+'\t'+key
            tmp_names.write(value+'\t'+key+'\n')
    run_metadata = False
    if 'y' in ARGS.metadata_run.lower() or 'y' in ARGS.roary_run.lower():
        run_metadata = True
        if 'y' in ARGS.roary_run.lower():
            print '\nRun roary requested, but we first need to run the'+\
                   ' metadata analysis.  Let\'s go!'
    if run_metadata:
       #summary_frames will store all of the metaDataFrames herein
        summary_frames = []
        n_isos = len(isos)

        #Kraken set at 2 threads, so 36 processes can run on 72 CPUs
        #Create a pool 'p' of size based on number of isolates (n_isos)
        if n_isos <= ARGS.threads//2:
            p = Pool(n_isos)
        else:
            p = Pool(ARGS.threads//2)
        print '\nRunning kraken on the assemblies (SPAdes contigs.fa files):'
        results_k_cntgs = p.map(kraken_contigs_multiprocessing, isos)
        #concat the dataframe objects
        res_k_cntgs = pd.concat(results_k_cntgs, axis=0)
        print '\nKraken_contigs results gathered from kraken on contigs...'

        #Multiprocessor retrieval of kraken results on reads.  Single thread
        #per job.
        if n_isos <= ARGS.threads:
            p = Pool(n_isos)
        else:
            p = Pool(ARGS.threads)
        results_k_reads = p.map(kraken_reads_multiprocessing, isos)
        #concat the dataframe objects
        res_k_reads = pd.concat(results_k_reads, axis=0)
        print 'Kraken_reads results gathered from kraken.tab files...'

        #Multiprocessor retrieval of contig metrics.  Single process
        #per job.
        results_metrics_contigs = p.map(metricsContigs_multiprocessing, isos)
        res_m_cntgs = pd.concat(results_metrics_contigs, axis=0)
        print 'Contig metrics gathered using \'fa -t\'...'

        #Multiprocessor retrieval of read metrics.  Single process
        #per job.
        results_metrics_reads = p.map(metricsReads_multiprocessing, isos)
        res_m_reads = pd.concat(results_metrics_reads, axis=0)
        print 'Read metrics gathered from yield.tab files...'

        #Multiprocessor retrieval of abricate results. Single process
        #per job.
        results_abricate = p.map(abricate_multiprocessing, isos)
        res_all_abricate = pd.concat(results_abricate, axis=0)
        res_all_abricate.fillna('.', inplace=True)
        print 'Resistome hits gathered from abricate.tab files...'

        #append the dfs to the summary list of dfs
        summary_frames.append(res_k_cntgs)
        summary_frames.append(res_k_reads)
        summary_frames.append(res_m_cntgs)
        summary_frames.append(res_m_reads)
        summary_frames.append(res_all_abricate)

        #These next steps build up the metadata not yet obtained
        #(via mulitprocesses above), also replace the dm-matrix short names
        #with original names

        #Let's store the metadata for each isolate in summary_isos
        summary_isos = []

        #Let's populate summary_isos above, isolate by isolate (in series)
        c = 0
        for iso in isos:
            iso_df = []
            if ARGS.excel_spreadsheet != None:
                #need to remove the suffixes so as to be able to match LIMS
                #excel metadata.
                if '-' in iso:
                    iso_nosuffix = iso.split('-')
                    iso_nosuffix = iso_nosuffix[0]+'-'+iso_nosuffix[1]
                else:
                    iso_nosuffix = iso
                submitter = xls_table.loc[iso_nosuffix, 'Submitter']
                maldi = xls_table.loc[iso_nosuffix,
                                      'Species identification (MALDI-TOF)']
                sp_id_subm = xls_table.loc[iso_nosuffix,
                                           'Species identification (Subm. lab)']
#                PCR_resgenes = xls_table.loc[iso_nosuffix,
#                                             'Final resist genes/alleles']
#                PCR_resgenes = str(PCR_resgenes)
#                PCR_resgenes = PCR_resgenes.replace(',', ';')
                lims = {'sp_LIMS_MALDI-Tof': maldi, 'sp_LIMS_SubmLab':
                        sp_id_subm, 'LIMS_Submitter': submitter}#,
                        #'Final_resgenes_PCR': PCR_resgenes}
                lims_df = pd.DataFrame([lims], index=[iso])
                while c < 1:
                    print 'LIMS metadata added to collection...'
                    c += 1
                iso_df.append(lims_df)
            sample = Isolate(iso)
            short_id = iso_ID_trans[iso]
            species_cntgs = res_k_cntgs.loc[iso, 'sp_krkn1_cntgs']
            species_reads = res_k_reads.loc[iso, 'sp_krkn1_reads']
            if species_cntgs == species_reads:
                species = species_cntgs
            else:
                species = 'indet'
            mlst_df = sample.mlst(species, sample.assembly())
            iso_df.append(mlst_df)
            species_consensus = {'sp_krkn_ReadAndContigConsensus':species}
            species_cons_df = pd.DataFrame([species_consensus], index=[iso])
            iso_df.append(species_cons_df)
            if new_ids != None:
                if iso in new_ids:
                    new_iso = {'0_new':'yes'}
                    new_iso_df = pd.DataFrame([new_iso], index=[iso])
                    iso_df.append(new_iso_df)
            iso_df_pd = pd.concat(iso_df, axis=1)
            summary_isos.append(iso_df_pd)
        print 'Remaining isolate data gathered (mlst, species consensus,'+\
              ' flag-if-new)...'
        #Glue the isolate by isolate metadata into a single df
        summary_isos_df = pd.concat(summary_isos)
        #Glue the dataframes built during multiprocessing processes
        summary_frames_df = pd.concat(summary_frames, axis=1)
        #Finish up with everything in one table!
        metadata_overall = pd.concat([summary_isos_df, summary_frames_df],
                                     axis=1)
        metadata_overall.fillna('.', inplace=True)
        print '\nMetadata super-matrix:'
        print metadata_overall
        #Write this supermatrix (metadata_overall) to csv and tab/tsv
        csv = base+'_metadataAll.csv'
        tsv = base+'_metadataAll.tab'
        json = base+'_metadataAll.json'
        metadata_overall.to_csv(csv, mode='w', index=True, index_label='name')
        metadata_overall.to_csv(tsv, mode='w', sep='\t', index=True,
                                index_label='name')
        metadata_overall.to_json(json)
        print '\nMetadata super-matrix for '+str(len(metadata_overall.index))+\
              ' isolates written to '+csv+' and '+tsv+'.'
        #Populate the isos_grouped_by_cons_spp dict with isolate IDs
        for k, v in zip(metadata_overall['sp_krkn_ReadAndContigConsensus'],
                        metadata_overall.index):
            isos_grouped_by_cons_spp[k.replace(' ', '_')].append(v)

    #Run andi?
    if 'y' in ARGS.andi_run.lower():
        #Run andi
        andi_mat = 'andi_'+ARGS.model_andi_distance+'dist_'+base+'.mat'
        andi_c = 'nice andi -j -m '+ARGS.model_andi_distance+' -t '+\
                  str(ARGS.threads)+' '+assembly_tempdir+'/*_contigs.fa > '+\
                  andi_mat
        print '\nRunning andi with: \''+andi_c+'\''
        os.system(andi_c)

        #Read in the andi dist matrix, convert to lower triangle
        dm = read_file_lines(andi_mat)[1:]
        dm = lower_tri(dm)
        #Correct the names in the matrix
        for iso in isos:
            #Could do it this way, but this is slower than a nested loop
            #dm.names[dm.names.index(iso_ID_trans[iso])] = iso
            #real	0m9.417s
            #user	1m18.576s
            #sys	0m2.620s
            #Nested loop is faster
            for i in range(0, len(dm.names)):
                #iso_ID_trans[iso] is the short_id
                if dm.names[i] == iso_ID_trans[iso]:
                    dm.names[i] = iso
            #real	0m8.789s
            #user	1m14.637s
            #sys	0m2.420s

        #From the distance matrix in dm, infer the NJ tree
        constructor = DistanceTreeConstructor()
        njtree = constructor.nj(dm)
        njtree.rooted = True
        Phylo.write(njtree, 'temp.tre', 'newick')
        t = Tree('temp.tre', format=1)
        #Get rid of negative branch lengths (an artefact, not an error, of NJ)
        for node in t.traverse():
            node.dist = abs(node.dist)
        t.set_outgroup(t.get_midpoint_outgroup())
        t_out = base+'_andi_NJ_'+ARGS.model_andi_distance+'dist.nwk.tre'
        t.write(format=1, outfile=t_out)
        print 'Final tree (midpoint-rooted, NJ under '+\
               ARGS.model_andi_distance+' distance) looks like this:'
        #Print the ascii tree
        print t
        #Remove the temp.tre
        os.remove('temp.tre')
        print 'Tree (NJ under '+ARGS.model_andi_distance+\
              ' distance, midpoint-rooted) written to '+t_out+'.'

    #Run roary?
    if 'y' in ARGS.roary_run.lower():
        roary_keepers = [
                        "accessory.header.embl",
                        "accessory.tab",
                        "accessory_binary_genes.fa",
                        "accessory_binary_genes.fa.newick",
                        "accessory_binary_genes_midpoint.nwk.tre",
                        "accessory_graph.dot",
                        "blast_identity_frequency.Rtab",
                        "clustered_proteins",
                        "core_accessory.header.embl",
                        "core_accessory.tab",
                        "core_accessory_graph.dot",
                        "core_gene_alignment.aln",
                        "gene_presence_absence.Ltab.csv",
                        "gene_presence_absence.Rtab",
                        "gene_presence_absence.csv",
                        "number_of_conserved_genes.Rtab",
                        "number_of_genes_in_pan_genome.Rtab",
                        "number_of_new_genes.Rtab",
                        "number_of_unique_genes.Rtab",
                        "pan_genome_reference.fa",
                        "pan_genome_sequences",
                        "summary_statistics.txt"
                        ]
        params = [(i, 'prokka') for i in isos if not
                  os.path.exists('prokka/'+i)]
        if len(params) > 0:
            print '\nRunning prokka:'
            if len(params) <= ARGS.threads//2:
                p = Pool(len(params))
            else:
                p = Pool(ARGS.threads//2)
            p.map(prokka, params)
        else:
            print '\nProkka files already exist. Let\'s move on to '+\
                  'the roary analysis...'

        #Run Roary on the species_consensus subsets.
        print 'Now, let\'s run roary!'
        for k, v in isos_grouped_by_cons_spp.items():
            print k, v
            n_isos = len(v)
            if n_isos > 1:
                shutil.rmtree(base+'_'+k+'_roary', ignore_errors=True)
                roary(base, k,
                      ' '.join(['prokka/'+iso+'/*.gff' for iso in v]))
                roary_genes = pd.read_table(base+'_'+k+
                                            '_roary/gene_presence_absence.' +\
                                            'Rtab',
                                            index_col=0, header=0)
                roary_genes = roary_genes.transpose()
                roary_genes.to_csv(base+'_'+k+
                                   '_roary/gene_presence_absence.Ltab.csv',
                                   mode='w', index=True, index_label='name')
                if n_isos > 2:
                    t = Tree(base+'_'+k+
                             '_roary/accessory_binary_genes.fa.newick',
                             format=1)
                    #Get rid of negative branch lengths (an artefact,
                    #not an error, of NJ)
                    for node in t.traverse():
                        node.dist = abs(node.dist)
                    t.set_outgroup(t.get_midpoint_outgroup())
                    t_out = base+'_'+k+\
                            '_roary/accessory_binary_genes_midpoint.nwk.tre'
                    t.write(format=1, outfile=t_out)
                    print '\nWritten midpoint-rooted roary tree.\n'
                    wd = os.getcwd()
                    os.chdir(base+'_'+k+'_roary')
                    for f_name in glob.glob('*'):
                        if f_name not in roary_keepers:
                            shutil.rmtree(f_name, ignore_errors=True)
                            os.remove(f_name)
                    os.chdir(wd)
                if n_isos <= 2:
                    print 'Need more than two isolates to have a meaningful '+\
                          'pangenome tree. No mid-point rooting of the ' +\
                          'pangenome tree performed.'
                wd = os.getcwd()
                os.chdir(base+'_'+k+'_roary')
                os.system('python ../collapseSites.py -f core_gene_alignment.aln -i fasta -t '+str(ARGS.threads))
                if os.path.exists('core_gene_alignment_collapsed.fasta'):
                    os.system('FastTree -nt -gtr < core_gene_alignment_collapsed.fasta > core_gene_FastTree_SNVs.tre')

                    #calc pairwise snp dist and write to file
                    with open('core_gene_alignment_collapsed.fasta', 'r') as inf:
                        aln = AlignIO.read(inf, 'fasta')
                        pairs = []
                        for i in range(0,len(aln)):
                            lst = [(aln, i, j) for j in range(0, i+1)]
                            pairs.append(lst)
                        if len(pairs) <= ARGS.threads:
                            p = Pool(len(pairs))
                        else:
                            p = Pool(ARGS.threads)
                        print 'Running pw comparisons in parallel...'
                        result = p.map(pw_calc, pairs)
                        summary = pd.concat(result, axis=0)
                        summary.fillna('', inplace=True)
                        with open('core_gene_alignment_SNV_distances.tab', 'w') as distmat:
                            summary.to_csv(distmat, mode='w', sep='\t', index=True, index_label='name')

                #convert roary output to fripan compatible
                os.system('python ../roary2fripan.py '+base+'_'+k)
                roary2fripan_strains_file = pd.read_table(base+'_'+k+
                                                          '.strains',
                                                          index_col=0,
                                                          header=0)
                info_list = []
                info_list.append(roary2fripan_strains_file)
                info_list.append(metadata_overall.loc[v, :])
                strains_info_out = pd.concat(info_list, axis=1)
                strains_info_out.to_csv(base+'_'+k+'.strains', mode='w',
                                        sep='\t', index=True,
                                        index_label='ID')
                print 'Updated '+base+'_'+k+'.strains with all metadata.'
                os.system('cp '+base+'_'+k+'* ~/public_html/fripan')
                os.chdir(wd)
            else:
                print 'Only one isolate in '+k+'. Need at least 2 isolates '+\
                      'to run roary.  Moving on...'

    #Delete the tempdirs created during the run
    if 'y' in ARGS.delete_tempdirs.lower():
        shutil.rmtree(assembly_tempdir, ignore_errors=True)
        print '\nDeleted tempdir '+assembly_tempdir+'.'
    else:
        print '\nTempdir '+assembly_tempdir+' not deleted.'

    print '\nRun finished.'
    print '\nExplore your results with phandango and FigTree.'
    print 'https://jameshadfield.github.io/phandango/'
    print 'http://tree.bio.ed.ac.uk/software/figtree/'
    print '\nContact the author at mark.schultz@unimelb.edu.au'
    print 'Thanks for using \''+VERSION+'\'.\n'


if __name__ == '__main__':
    main()
