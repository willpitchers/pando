#!/usr/bin/env python3

import os

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
        cmd_kraken = 'nice kraken --threads 2 --db /home/linuxbrew/db/kraken//microbekraken'+\
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