![pando](http://upliftconnect.com/wp-content/uploads/2016/03/pando-trees-1.jpg)
# pando
With a file of isolate IDs as input, compile all results (abricate, kraken, mlst, contig and read metrics, LIMS metadata) into a single super-matrix.  Options (all default to 'don't run'): infer an NJ tree using the alignment-free Andi phylogenomic software; gather metdata for all isolates; run roary on the isolate set (have a prokka folder in the rundir to save computation time).  Use a method of choice for displaying the metadata next to the tree (e.g. we recommend [phandango](https://jameshadfield.github.io/phandango/) (web based) and [FigTree](http://tree.bio.ed.ac.uk/software/figtree/) (installed locally on your machine)).  Can Use it to run roary too.  Can switch on/off roary or andi or both.  

### Example run command on MDU servers:
`time nice python pando.py -i file_of_mdu-isolate_IDs.txt -n 2016-15949 2016-15442 -x CPE_ongoing_20160801.xlsx -m y`

### Minimum run command:
`python pando.py -i file_of_mdu-isolate_IDs.txt`

### Get help:
```
python pando.py -h
usage: pando.py [-h] -i MDU_READ_IDS [-n NEW_IDS [NEW_IDS ...]] [-w WGS_QC]
                [-d DELETE_TEMPDIRS] [-t THREADS] [-a ANDI_RUN] [-r ROARY_RUN]
                [-m METADATA_RUN] [-s MODEL_ANDI_DISTANCE] [-c PERCENT_CUTOFF]
                [-x EXCEL_SPREADSHEET]

Run exploratory analyses.

optional arguments:
  -h, --help            show this help message and exit
  -i MDU_READ_IDS, --mdu_read_IDs MDU_READ_IDS
                        One MDU-ID per line (required). Put in same folder as
                        run folder.
  -n NEW_IDS [NEW_IDS ...], --new_IDs NEW_IDS [NEW_IDS ...]
                        Enter IDs (space delimited) that you wish to 'flag-if-
                        new' in the final table.
  -w WGS_QC, --wgs_qc WGS_QC
                        Path to WGS QC. Default='/mnt/seq/MDU/QC/'
  -d DELETE_TEMPDIRS, --delete_tempdirs DELETE_TEMPDIRS
                        Delete tempdirs created during run? Default='yes'.
  -t THREADS, --threads THREADS
                        Number of threads, default='72'
  -a ANDI_RUN, --andi_run ANDI_RUN
                        Run andi phylogenomic analysis? Default='no'
  -r ROARY_RUN, --roary_run ROARY_RUN
                        Run roary pangenome analysis? Default='no'
  -m METADATA_RUN, --metadata_run METADATA_RUN
                        Gather metadata for all isolates? Default='no'
  -s MODEL_ANDI_DISTANCE, --model_andi_distance MODEL_ANDI_DISTANCE
                        Substitution model. 'Raw', 'JC', or 'Kimura'.
                        Default='JC'.
  -c PERCENT_CUTOFF, --percent_cutoff PERCENT_CUTOFF
                        For abricate, call the gene 'present' if greater than
                        this value and 'maybe' if less than this value.
                        Default=95. NB: 100 percent= '100', not '1'.
  -x EXCEL_SPREADSHEET, --excel_spreadsheet EXCEL_SPREADSHEET
                        Parse excel spreadsheet of metadata (.xlsx format)to
                        extract LIMS data. The data must start on line 5
                        (1-based indexing, as per Excel line numbers), which
                        is default for LIMS, and contain columns with the
                        labels 'MDU sample ID', 'Species identification
                        (MALDI-TOF)', 'Species identification (Subm. lab)' and
                        'Submitter'.
```

