![pando](http://upliftconnect.com/wp-content/uploads/2016/03/pando-trees-1.jpg)
# pando
With a file of isolate IDs as input, compile all results (abricate, kraken, mlst) into a single table and infer an NJ tree using the alignment-free Andi phylogenomic software.  Receive the tree and results table as an email at the end of the run.

### Example run command:
`time nice python pando.py -i isolates_fileOfFileNames.txt -e user@domain.edu.au -j jobNumberX -d y`

### Get help:
```
python pando.py -h
usage: pando.py [-h] -i MDU_READ_IDS [-w WGS_QC] [-d DELETE_TEMPDIRS]
                [-t THREADS] [-a ANDI_RUN] [-m MODEL_ANDI_DISTANCE]
                [-c PERCENT_CUTOFF] -e EMAIL_ADDRESSES [EMAIL_ADDRESSES ...]
                -j JOB_NUMBER

CPE ongoing

optional arguments:
  -h, --help            show this help message and exit
  -i MDU_READ_IDS, --mdu_read_ids MDU_READ_IDS
                        One MDU-ID per line. Put in same folder as run folder.
  -w WGS_QC, --wgs_qc WGS_QC
                        Path to WGS QC. Default '/mnt/seq/MDU/QC/'
  -d DELETE_TEMPDIRS, --delete_tempdirs DELETE_TEMPDIRS
                        Delete tempdirs created during run? Default = 'no'.
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
```

