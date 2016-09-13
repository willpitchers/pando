# Pando, analysis workflow
In its current state, pando.py relies on the MDU file structure and initial QC results.  It also depends on a number of softwares installed on the mdu-servers.  To use it outside of the MDU file system, some parts of the code will need to be re-written.  

## 1. Intitiating a job.

A LIMS request will be submitted to Bioinformatics as an excel spreadsheet, for CPE with a file name CPE_ongoing_YYYMMDD.xlsx (YYYYMMDD = Year, Month, Date).  The file contains a set of CPE isolates in a two year rolling window.  Open the file and look for the red highlighted rows – these are the new isolates for this request.

An email will also be sent from Molecular to Bioinformatics stating that the reads are now available for analysis.  Bioinformatics must perform initial QC on the reads before pando can be run successfully.  They will let the analyst know when this is in progress or has been completed.  

Copy the excel spreadsheet to an existing folder on the server:
```
scp /Users/local_username/Desktop/MDU_DAMG/MDU/CPE_ongoing/CPE_ongoing_2010913.xlsx remote_username@zaphod.mdu.unimelb.edu.au:/home/remote_username/jobs/mdu/CPE_ongoing/pando
```

Download the required pando scripts from github (the latter two scripts are only required if you wish to perform a roary analysis [see section 2 below]; thanks to Jason Kwong for the `roary2fripan.py` script):
```
wget raw.githubusercontent.com/MDU-PHL/pando/master/pando.py
wget raw.githubusercontent.com/MDU-PHL/pando/master/collapseSites.py
wget raw.githubusercontent.com/kwongj/roary2fripan/master/roary2fripan.py
```

Test it with:
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


Copy the first column of MDU-IDs from the spreadsheet and paste into a single column text file, with no header (named isolates_YYYMMDD.txt).

```
nano isolates_2010913.txt
```

Start a screen session to run the job, named using `-S`:
```
screen -S pando
```

Start a screen log to save the stdout to text file using the keystroke
Control+a, shift+H.  Toggle the screenlog.0 off with the same keystroke.   

## 2. Running the job request using pando
Test the setup using:
```
time nice python pando.py -i isolates_20160913.txt -n 2016-18848 2016-18844 2016-18808 2016-18697 2016-18648 -x CPE_ongoing_2010913.xlsx
```

There a three options to run pando:

<i>i</i>) `-m y` a metadata only run (will create a metadata super-matrix for all isolates in the analysis and report those that were missing). This step is fast.<br>
<i>ii</i>) `-a y` an andi run (to get the NJ tree). This step is slow. <br>
<i>iii</i>) `-r y` a roary run to get the pangenome and associated fripan files (do a `mkdir /home/username/public_html/fripan` before running this option, and clone the fripan code from https://github.com/drpowell/FriPan into this folder).  This step is very slow.<br>  
Assess the load on the server using `htop` and choose the number of threads acccordingly (`72` by default)<br>

To do a full run:
```
time nice python pando.py -i isolates_20160913.txt -n 2016-18848 2016-18844 2016-18808 2016-18697 2016-18648 -x CPE_ongoing_2010913.xlsx -m y -a y -r y
```
Whilst running, exit the screen using control+a, d.  Check the progress by `tail screenlog.0`, `htop -u username` or going back into the screen with `screen -r pando`. <br>
Kill the screen process ID (`PID`) using `top` or enter the screen, hit control+c then control+d (this latter option might not work if `pando` is in the middle of a multiprocessing loop).  Rename the `screenlog.0` file at the end of the run if you would like to store it for future reference (else it may be overwritten). 
