![pando](http://upliftconnect.com/wp-content/uploads/2016/03/pando-trees-1.jpg)
# pando
With a file of isolate IDs as input, compile all results (abricate, kraken, mlst, contig and read metrics, LIMS metadata) into a single super-matrix.  Infer an NJ tree using the alignment-free Andi phylogenomic software.  Receive the tree and results table as an email at the end of the run.  Use a method of choice for displaying the metadata next to the tree (e.g. we recommend [phandango](https://jameshadfield.github.io/phandango/) (web based) and FigTree (installed locally on your machine))

### Example run command:
` time nice python pando.py -i isos.txt -e recipient@unimelb.edu.au -j job_2467 -n 2016-15949 2016-15442 -x CPE_ongoing_20160801.xlsx`

### Get help:
```
python pando.py -h
usage: pando.py [-h] -i MDU_READ_IDS [-n NEW_IDS [NEW_IDS ...]] [-w WGS_QC]
                [-d DELETE_TEMPDIRS] [-t THREADS] [-a ANDI_RUN]
                [-m MODEL_ANDI_DISTANCE] [-c PERCENT_CUTOFF] -e
                EMAIL_ADDRESSES [EMAIL_ADDRESSES ...] -j JOB_NUMBER -x
                EXCEL_METADATA

Rrun expoloratory analyses

optional arguments:
  -h, --help            show this help message and exit
  -i MDU_READ_IDS, --mdu_read_IDs MDU_READ_IDS
                        One MDU-ID per line. Put in same folder as run folder.
  -n NEW_IDS [NEW_IDS ...], --new_IDs NEW_IDS [NEW_IDS ...]
                        Enter IDs (space delimited) that you wish to have
                        flagged as 'new' in the final table.
  -w WGS_QC, --wgs_qc WGS_QC
                        Path to WGS QC. Default '/mnt/seq/MDU/QC/'
  -d DELETE_TEMPDIRS, --delete_tempdirs DELETE_TEMPDIRS
                        Delete tempdirs created during run? Default = 'yes'.
  -t THREADS, --threads THREADS
                        Number of threads, default='72'
  -a ANDI_RUN, --andi_run ANDI_RUN
                        Run andi phylogenomic analysis? Default='yes'
  -m MODEL_ANDI_DISTANCE, --model_andi_distance MODEL_ANDI_DISTANCE
                        Substitution model. 'Raw', 'JC', or 'Kimura'. Default
                        = 'JC'.
  -c PERCENT_CUTOFF, --percent_cutoff PERCENT_CUTOFF
                        For abricate, call the gene 'present' if greater than
                        this value and 'maybe' if less than this value.
                        Default = 95. NB: 100 percent = '100', not '1'.
  -e EMAIL_ADDRESSES [EMAIL_ADDRESSES ...], --email_addresses EMAIL_ADDRESSES [EMAIL_ADDRESSES ...]
                        Email addresses to send results (comma or space
                        delimited)
  -j JOB_NUMBER, --job_number JOB_NUMBER
                        Enter the MDU job number (no spaces).
  -x EXCEL_METADATA, --excel_metadata EXCEL_METADATA
                        Parse excel spreadsheet of metadata.
```

