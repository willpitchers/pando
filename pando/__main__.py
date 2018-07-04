#!/usr/bin/env python3


'''

Uses python3.
Email: dr.mark.schultz@gmail.com
Github: schultzm
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.


Read in metadata from Excel spreadsheet.
Return a supermatrix of QC metadata and a tree.
Specific for MDU folder structures and QC
Optionally run roary analysis or tree
Email: dr.mark.schultz@gmail.com
Github: https://github.com/schultzm

Acknowledgements:
Will Pitchers (usage, bug finding, feature requests)
Torsten Seemann (pando relies heavily on Torsten's tools)
Anders Goncalves da Silva (advice on multiprocessing)
Dieter Bulach (advice on setting up fripan on the MDU server)
Jason Kwong (roary2fripan, advised to run Kraken on the contigs)
Susan Ballard (ongoing feature requests)
Tim Stinear (advised to use more-than-kraken to bin isolates into species)
Andrew Page (for all roary tools)
Ben Howden (feature requests)

'''


import os
import argparse
import sys
import random
from subprocess import Popen, PIPE
import shlex
import shutil
import glob
from multiprocessing import Pool
import pandas as pd
from multiprocessing import cpu_count


# set up the arguments parser to deal with the command line input
PARSER = argparse.ArgumentParser(description='Run QC summary analysis.',
                                 prog='pando',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
SUBPARSER_ARGS1 = argparse.ArgumentParser(add_help=False)
SUBPARSER_ARGS1.add_argument('LIMS_request_sheet', help="Excel spreadsheet LIMS request.")
SUBPARSER_ARGS1.add_argument('-w', '--wgs_qc', help='Path to WGS QC',
                    default='/mnt/seq/MDU/QC/', required=False)
SUBPARSER_ARGS1.add_argument('-N', '--Nullarbor_folders', help='Run this on Nullarbor outdir?',
                            default=False, action="store_true", required=False)
SUBPARSER_ARGS1.add_argument('-k', '--keep_tempdirs', help='Keep tempdirs created\
                    during run?', default=False, action='store_true',
                    required=False)
SUBPARSER_ARGS1.add_argument('-d', '--kraken_db', help='Path to Kraken db',
                             default=os.environ['KRAKEN_DEFAULT_DB'])
SUBPARSER_ARGS1.add_argument("-t", "--threads", help='Number of threads',
                    default=cpu_count(), type=int, required=False)
SUBPARSER_ARGS1.add_argument('-a', '--andi_run', help='Run andi phylogenomic analysis?\
                    ', default=False, action='store_true', required=False)
SUBPARSER_ARGS1.add_argument('-r', '--roary_run', help='Run roary pangenome analysis?\
                    ', default=False, action='store_true', required=False)
SUBPARSER_ARGS1.add_argument('-m', '--metadata_run', help='Gather metadata for all\
                    isolates?', default=True, action='store_false')
SUBPARSER_ARGS1.add_argument('-s', '--model_andi_distance', help='Substitution model.\
                    \'Raw\', \'JC\', or \'Kimura\'.',
                    default='JC', required=False)
SUBPARSER_ARGS1.add_argument('-c', '--percent_cutoff', help='For abricate, call the\
                    gene \'present\' if greater than this value and \'maybe\'\
                    if less than this value.', default=95, type=int, required=False)
SUBPARSER_MODULES = PARSER.add_subparsers(title="Sub-commands help",
                                          help="",
                                          metavar="",
                                          dest="subparser_name")
SUBPARSER_MODULES.add_parser("run",
                             help="""Run MDU AMR request.""",
                             description="Start the run.",
                             parents=[SUBPARSER_ARGS1],
                             formatter_class=argparse.ArgumentDefaultsHelpFormatter)
SUBPARSER_MODULES.add_parser("version", help="Print version.",
                             description="Print version.")

ARGS = PARSER.parse_args()


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
        isolate_qc_contigs = os.path.join(ARGS.wgs_qc, self.ID, 'contigs.fa')
        if os.path.exists(isolate_qc_contigs):
            return isolate_qc_contigs
        else:
            print(isolate_qc_contigs+' does not exist')

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
        metrics = dict(list(zip(metrics[0], metrics[1])))
        os.system('rm '+self.ID+'_metrics.txt')
        metrics_df = pd.DataFrame([metrics], index=[self.ID])
        return metrics_df

    def get_yield(self):
        '''
        Return the read metrics stored in yield.tab QC file.
        '''
        yield_data = [line.rstrip().split('\t') for line in open(ARGS.wgs_qc+\
                      self.ID+'/'+YIELD_FILE).readlines()]
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
        Get abricate results or run abricate.
        '''
        wd = ARGS.wgs_qc+self.ID
        abricate_path = wd+'/abricate.tab'
        if os.path.exists(abricate_path):
            ab_data = pd.read_table(abricate_path, sep='\t', header=0)
        else: #run abricate.
            abricate_outfolder = 'abricate/'+self.ID
            abricate_outfile = abricate_outfolder+'/abricate.tab'
            contigs_path = wd+'/contigs.fa'
            if os.path.exists(abricate_outfile):
                pass
            else:
                if os.path.exists(contigs_path):
                    print('running abricate for', self.ID)
                    os.system('mkdir -p '+abricate_outfolder)
                    os.system('abricate '+contigs_path+' > '+abricate_outfile)
            ab_data = pd.read_table(abricate_outfile, sep='\t', header=0)
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
        import itertools
        ab_results = dict(itertools.chain(iter(y.items()), iter(m.items())))
        #Convert to pandas dataframe
        abricate_results = pd.DataFrame([ab_results], index=[self.ID])
        return abricate_results

    def mlst(self, species, assembly):
        '''
        Store the MLST. If the species is listed in FORCE_MLST_SCHEME redo
        MLST.  If the mlst.tab exists and the species is not in FORCE, use
        the existing mlst.tab if running in /mnt/seq/MDU/QC or mlst2.tab
        if using a Nullarbor folder structure.
        '''
        if species in FORCE_MLST_SCHEME:
            cmd = 'mlst --legacy --scheme '+FORCE_MLST_SCHEME[species]+' --quiet '+\
                   assembly
            args_mlst = shlex.split(cmd)
            #should write a proc function to do this call as it's used often
            proc = Popen(args_mlst, stdout=PIPE)
            output = proc.stdout.read().decode('UTF-8')
            mlst = output.rstrip().split('\n')
            mlst = [line.strip().split('\t') for line in [_f for _f in mlst if _f]]
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
            mlst_tab = ARGS.wgs_qc+self.ID+'/'+MLST_FILE
            run_mlst_again = False
            if os.path.exists(mlst_tab):
                mlst = [line.rstrip().split('\t') for line in
                        open(mlst_tab).readlines()][0]
                if 'SCHEME' in mlst:
                    #This would mean the old version of MLST was used
                    run_mlst_again = True
                mlst_formatted_dict = {'MLST_Scheme': mlst[1],
                                       'MLST_ST': mlst[2]}
                k = 1
                for i in range(3, len(mlst)):
                    mlst_formatted_dict['MLST_Locus'+str(k)] = mlst[i]
                    k += 1
            if run_mlst_again or os.path.exists(mlst_tab) == False:
                cmd = 'mlst --quiet '+assembly
                args_mlst = shlex.split(cmd)
                proc = Popen(args_mlst, stdout=PIPE)
                output = proc.stdout.read().decode('UTF-8')
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
        cmd_sort = 'sort -k 1 -g -r'
        cmd_head = 'head -3'
        #Split the cmds using shlex, store in args
        args_grep = shlex.split(cmd_grep)
        args_sort = shlex.split(cmd_sort)
        args_head = shlex.split(cmd_head)

        #Pipe the output of one args to another
        proc1 = Popen(args_grep, stdout=PIPE)
        proc2 = Popen(args_sort, stdin=proc1.stdout, stdout=PIPE)
        proc3 = Popen(args_head, stdin=proc2.stdout, stdout=PIPE)

        output = proc3.stdout.read().decode('UTF-8')
        kraken = output.rstrip().split('\n')
        kraken = [line.strip().split('\t') for line in [_f for _f in kraken if _f]]
        return kraken

    def kraken_contigs(self):
        '''
        Get the kraken best hit from assemblies.
        '''
        #Pipe these commands together
        cmd_kraken = 'nice kraken --threads 2 --db '+\
                     os.path.abspath(ARGS.kraken_db)+\
                     ' --fasta-input '+ARGS.wgs_qc+'/'+self.ID+'/contigs.fa'
        cmd_krk_r = 'kraken-report'
        cmd_grep = "grep -P '\tS\t'"
        cmd_sort = 'sort -k 1 -g -r'
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

        output = proc5.stdout.read().decode('UTF-8')
        kraken = output.rstrip().split('\n')
        kraken = [line.strip().split('\t') for line in [_f for _f in kraken if _f]]
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
    import string
    chars = string.ascii_uppercase + string.digits
    return ''.join(random.SystemRandom().choice(chars) for _ in range(10))

def make_tempdir():
    '''
    Make a temporary directory.
    '''
    import tempfile
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
        lower_triangle.append(list(map(float, i[1:k])))
        names.append(i[0])
        k += 1
    from Bio.Phylo.TreeConstruction import _DistanceMatrix
    matrix = _DistanceMatrix(names, lower_triangle)
    return matrix

def get_isolate_request_IDs(ID_file):
    '''
    Reads in the MDU IDs from the request IDs file and returns IDs as a list.
    ID file must contain only one ID per line.
    '''
    df = pd.read_excel(ID_file, skiprows=0, index_col=0)
    return df

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
    if len(avail) > 0:
        print('\nFound folders with the following IDs:')
        print('\n'.join(avail))
    not_avail = sorted(list(set(not_avail)))
    if len(not_avail) > 0:
        outf = os.path.splitext(ARGS.LIMS_request_sheet)[0]+'_not_found.txt'
        with open(outf, 'w') as not_found:
            not_found.write('\n'.join(not_avail))
        print('\nCould not find IDs:\n', '\n'.join(not_avail)+'\n')
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

def prokka(params):
    '''
    Run prokka on the isolate. Unpack params to get iso and assembly_tempdir.
    '''
    (x, y) = params
    ID = Isolate(x)
    cmd = ID.prokka_contigs() % y
    print(cmd)
    os.system(cmd)

def roary(base, sp, gffs):
    ''''
    Run roary on the gff files output by prokka.
    '''
    cmd = 'nice roary -f '+base+'_'+sp+'_roary -e -n -v -z -p '+\
          str(ARGS.threads)+' '+gffs
    print('\nRunning roary for '+sp+' with the command:\n'+cmd)
    os.system(cmd)

def pw_calc(aln_seq_coords):
    '''
    Calculate the absolute distance between two sequences, in a range of seqs.
    '''
    df_dist = []
    for m in aln_seq_coords:
        (aln, i, j) = m
        print(i, j)
        name = aln[i].id
        x = len([y for y in zip(aln[i].seq,
                aln[j].seq) if len(set(y)) > 1 and 'N' not in set(y) and '-' not in set(y) and '?' not in set(y)])
        df = pd.DataFrame([{aln[j].id: x}], index=[name])
        df_dist.append(df)
    return pd.concat(df_dist, axis=1, sort=False)


def main():
    global YIELD_FILE
    global MLST_FILE
    global FORCE_MLST_SCHEME
    #Set up the file names for Nullarbor folder structure
    YIELD_FILE = 'yield.tab'
    MLST_FILE = 'mlst.tab'


    #Add MLST schemes to force their usage if that species is encountered
    #Only force schemes if there are two (e.g., A baumannii and E coli)
    FORCE_MLST_SCHEME = {"Acinetobacter baumannii": "abaumannii_2",
                         "Campylobacter jejuni": "campylobacter",
                         #"Citrobacter freundii": "cfreundii",
                         #"Cronobacter": "cronobacter",
                         "Enterobacter cloacae": "ecloacae",
                         "Escherichia coli": "ecoli",
                         #"Klebsiella oxytoca": "koxytoca",
                         #"Klebsiella pneumoniae": "kpneumoniae",
                         #"Pseudomonas aeruginosa": "paeruginosa"
                         "Shigella sonnei": "ecoli",
                         "Salmonella enterica": "senterica",
                         "Vibrio cholerae": "vcholerae"
                        }


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
    if not ARGS.subparser_name:
        PARSER.print_help()
        sys.exit()


    elif ARGS.subparser_name == 'version':
        from .utils.version import Version
        Version()
        sys.exit()

    else:# ARGS.subparser_name == "run":
        if ARGS.Nullarbor_folders:
            print('Nullarbor folder structure selected.')
            YIELD_FILE = 'yield.clean.tab'
            MLST_FILE = 'mlst2.tab'

        EXCEL_OUT = (f"{os.path.splitext(os.path.basename(ARGS.LIMS_request_sheet))[0]}" \
                     f"_results.xlsx")

        if ARGS.threads > cpu_count():
            sys.exit(f'Number of requested threads must be less than {cpu_count()}.')

        print(str(ARGS.threads) +' CPU processors requested.')


        #Check if final slash in manually specified wgs_qc path
        if ARGS.wgs_qc[-1] != '/':
            print('\n-wgs_qc path is entered as '+ARGS.wgs_qc)
            print('You are missing a final \'/\' on this path.')
            print('Exiting now.\n')
            sys.exit()



        #i) read in the IDs from file
        xls_table = get_isolate_request_IDs(ARGS.LIMS_request_sheet)
        IDs = list(set(xls_table.index.values))

        #base should be a global, given that it is used in other functions too.
        base = os.path.splitext(ARGS.LIMS_request_sheet)[0]

        #ii) Return a folder path to the QC data for each available ID
        #    using a wildcard search of the ID in IDs in ARGS.wgs_qc path.
        iso_paths = isolates_available(IDs)
        #Drop the path and keep the folder name
        isos = [i.split('/')[-1] for i in iso_paths]

        #iii) make tempdir to store the temp_contigs there for 'andi' analysis.
        assembly_tempdir = make_tempdir()

        #vi) Copy contigs to become temp_contigs into tempdir, only if andi
        #requested.
        #Translation dict to store {random 9-character filename: original filename}
        iso_ID_trans = {}
        #Dict to store each isolate under each consensus species#####maybe delete
        from collections import defaultdict
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
            if ARGS.andi_run:
                cmd = 'ln -s '+assembly_path+' '+assembly_tempdir+'/'+short_id+\
                      '_contigs.fa'
                os.system(cmd)
                print('Creating symlink:', cmd)
        if len(list(iso_ID_trans.items())) > 0:
            with open(base+'_temp_names.txt', 'w') as tmp_names:
                print('\nTranslated isolate IDs:\nShort\tOriginal')
                for key, value in list(iso_ID_trans.items()):
                    print(value+'\t'+key)
                    tmp_names.write(value+'\t'+key+'\n')
        if ARGS.metadata_run:
           #summary_frames will store all of the metaDataFrames herein
            summary_frames = []
            n_isos = len(isos)
            if n_isos == 0:
                print('\nNo isolates detected in the path '+ARGS.wgs_qc+'.')
                print('Exiting now.\n')
                sys.exit()
            #Kraken set at 2 threads, so 36 processes can run on 72 CPUs
            #Create a pool 'p' of size based on number of isolates (n_isos)
            if n_isos <= ARGS.threads//2:
                p = Pool(n_isos)
            else:
                p = Pool(ARGS.threads//2)
            print('\nRunning kraken on the assemblies (SPAdes contigs.fa files):')
            results_k_cntgs = p.map(kraken_contigs_multiprocessing, isos)
            #concat the dataframe objects
            res_k_cntgs = pd.concat(results_k_cntgs, axis=0, sort=False)
            print('\nKraken_contigs results gathered from kraken on contigs...')

            #Multiprocessor retrieval of kraken results on reads.  Single thread
            #per job.
            if n_isos <= ARGS.threads:
                p = Pool(n_isos)
            else:
                p = Pool(ARGS.threads)
            results_k_reads = p.map(kraken_reads_multiprocessing, isos)
            #concat the dataframe objects
            res_k_reads = pd.concat(results_k_reads, axis=0)
            print('Kraken_reads results gathered from kraken.tab files...')

            #Multiprocessor retrieval of contig metrics.  Single process
            #per job.
            results_metrics_contigs = p.map(metricsContigs_multiprocessing, isos)
            res_m_cntgs = pd.concat(results_metrics_contigs, axis=0)
            print('Contig metrics gathered using \'fa -t\'...')

            #Multiprocessor retrieval of read metrics.  Single process
            #per job.
            results_metrics_reads = p.map(metricsReads_multiprocessing, isos)
            res_m_reads = pd.concat(results_metrics_reads, axis=0)
            print('Read metrics gathered from '+YIELD_FILE+' files...')

            #Multiprocessor retrieval of abricate results. Single process
            #per job.
            results_abricate = p.map(abricate_multiprocessing, isos)
            res_all_abricate = pd.concat(results_abricate, axis=0, sort=False)
            res_all_abricate.fillna('', inplace=True)
            print('Resistome hits gathered from abricate.tab files...')

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
                iso_df_pd = pd.concat(iso_df, axis=1)
                summary_isos.append(iso_df_pd)

            #Glue the isolate by isolate metadata into a single df
            summary_isos_df = pd.concat(summary_isos)
            #Glue the dataframes built during multiprocessing processes
            summary_frames_df = pd.concat(summary_frames, axis=1)
            #Finish up with everything in one table!
            metadata_overall = pd.concat([xls_table, summary_isos_df, summary_frames_df],
                                         axis=1, sort=False)

            metadata_overall.fillna('', inplace=True)
            metadata_overall.index.name = 'ISOLATE'
            print('\nMetadata super-matrix:')
            #Write this supermatrix (metadata_overall) to csv and tab/tsv
            csv = os.path.abspath(base+'_metadataAll.csv')
            tsv = os.path.abspath(base+'_metadataAll.tab')
            json = os.path.abspath(base+'_metadataAll.json')
            metadata_overall.to_csv(sys.stdout)
            writer = pd.ExcelWriter(EXCEL_OUT)
            metadata_overall.to_excel(writer,'Sheet 1', freeze_panes=(1, 1))
            writer.save()
            print(f"\nResults written to {os.path.abspath(EXCEL_OUT)}")

            for k, v in zip(metadata_overall['sp_krkn_ReadAndContigConsensus'],
                            metadata_overall.index):
                isos_grouped_by_cons_spp[k.replace(' ', '_')].append(v)

        #Run andi?
        if ARGS.andi_run:
            #Run andi
            andi_mat = 'andi_'+ARGS.model_andi_distance+'dist_'+base+'.mat'
            andi_c = 'nice andi -j -m '+ARGS.model_andi_distance+' -t '+\
                      str(ARGS.threads)+' '+assembly_tempdir+'/*_contigs.fa > '+\
                      andi_mat
            print('\nRunning andi with: \''+andi_c+'\'')
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
            from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
            constructor = DistanceTreeConstructor()
            njtree = constructor.nj(dm)
            njtree.rooted = True
            from Bio import Phylo
            Phylo.write(njtree, 'temp.tre', 'newick')
            from ete3 import Tree
            t = Tree('temp.tre', format=1)
            #Get rid of negative branch lengths (an artefact, not an error, of NJ)
            for node in t.traverse():
                node.dist = abs(node.dist)
            t.set_outgroup(t.get_midpoint_outgroup())
            t_out = base+'_andi_NJ_'+ARGS.model_andi_distance+'dist.nwk.tre'
            t.write(format=1, outfile=t_out)
            print('Final tree (midpoint-rooted, NJ under '+\
                   ARGS.model_andi_distance+' distance) looks like this:')
            #Print the ascii tree
            print(t)
            #Remove the temp.tre
            os.remove('temp.tre')
            print('Tree (NJ under '+ARGS.model_andi_distance+\
                  ' distance, midpoint-rooted) written to '+t_out+'.')

        #Run roary?
        if ARGS.roary_run:
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
                print('\nRunning prokka:')
                if len(params) <= ARGS.threads//2:
                    p = Pool(len(params))
                else:
                    p = Pool(ARGS.threads//2)
                p.map(prokka, params)
            else:
                print('\nProkka files already exist. Let\'s move on to '+\
                      'the roary analysis...')

            #Run Roary on the species_consensus subsets.
            print('Now, let\'s run roary!')
            for k, v in list(isos_grouped_by_cons_spp.items()):
                print(k, v)
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
                        from ete3 import Tree
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
                        print('\nWritten midpoint-rooted roary tree.\n')
                        wd = os.getcwd()
                        os.chdir(base+'_'+k+'_roary')
                        for f_name in glob.glob('*'):
                            if f_name not in roary_keepers:
                                shutil.rmtree(f_name, ignore_errors=True)
                                os.remove(f_name)
                        os.chdir(wd)
                    if n_isos <= 2:
                        print('Need more than two isolates to have a meaningful '+\
                              'pangenome tree. No mid-point rooting of the ' +\
                              'pangenome tree performed.')
                    wd = os.getcwd()
                    os.chdir(base+'_'+k+'_roary')
                    os.system('python ../collapseSites.py -f core_gene_alignment.aln -i fasta -t '+str(ARGS.threads))
                    if os.path.exists('core_gene_alignment_collapsed.fasta'):
                        os.system('FastTree -nt -gtr < core_gene_alignment_collapsed.fasta > core_gene_FastTree_SNVs.tre')

                        #calc pairwise snp dist and write to file
                        with open('core_gene_alignment_collapsed.fasta', 'r') as inf:
                            from Bio import AlignIO
                            aln = AlignIO.read(inf, 'fasta')
                            pairs = []
                            for i in range(0,len(aln)):
                                lst = [(aln, i, j) for j in range(0, i+1)]
                                pairs.append(lst)
                            if len(pairs) <= ARGS.threads:
                                p = Pool(len(pairs))
                            else:
                                p = Pool(ARGS.threads)
                            print('Running pw comparisons in parallel...')
                            result = p.map(pw_calc, pairs)
                            summary = pd.concat(result, axis=0, sort=False)
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
                    strains_info_out = pd.concat(info_list, axis=1, sort=False)
                    strains_info_out.to_csv(base+'_'+k+'.strains', mode='w',
                                            sep='\t', index=True,
                                            index_label='ID')
                    print('Updated '+base+'_'+k+'.strains with all metadata.')
                    os.system('cp '+base+'_'+k+'* ~/public_html/fripan')
                    os.chdir(wd)
                else:
                    print('Only one isolate in '+k+'. Need at least 2 isolates '+\
                          'to run roary.  Moving on...')

        #Keep the tempdirs created during the run
        if not ARGS.keep_tempdirs:
            shutil.rmtree(assembly_tempdir, ignore_errors=True)
            print('\nDeleted tempdir '+assembly_tempdir+'.')
        else:
            print('\nTempdir '+assembly_tempdir+' not deleted.')

        print('\nRun finished.')


