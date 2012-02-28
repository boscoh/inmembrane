# inmembrane

_inmembrane_ is a pipeline for proteome annotation to predict if a protein is exposed on the surface of a bacteria. 

As _inmembrane_ is a glue script that runs a lot of external dependencies, we have tried to make the parsing code for each program as concise as possible. The program is written to use as much of the standard Python library as possible, and specifically uses native data structures which provide an incredible amount of flexibility.

It is a fact of life for bioinformatics that new versions of basic tools changes output formats and API. We believe that it is an essential skill to rewrite parsers to handle the subtle but significant changes in different versions. Hopefully _inmembrane_ makes this as simple as possible.


## Execution

_inmembrane_ was written in Python 2.7. It takes a FASTA input file and runs  a number of external bioinformatic analysis programs on the sequences. It then collects the output to make the final analysis. It can be run in two modes. 

It can be run as a command-line program in the _inmembrane_ directory:  
     
    python inmembrane.py your_fasta_file

The other way of running imembrane.py is with a custom script, such as  `run_example.py`. If you open it with a **text editor**, you will see that the input file and the output file are filled in. The output file will be created for you. You can run it on the command-line: 

    python run_example.py

Or easier still, on most systems with Python installed, you simply double-click it in a file-manager.

## Output format

Given that this is a glue script that runs a lot of programs, there will be many places for error, and as such, to aid in debugging, we print out all the commands run in the script. These lines will be preceeded by the hash # character so that it can easily be commented out.

The output is CSV compatible. You can open it with EXCEL. It is also designed to be whitespace delineated. For instance, if no feature of interest is found for a sequence, then a '.' is put in the third column. The full name, which has variable number of tokens is put at the end. Here's an example:

    SPy_0008 , CYTOPLASM    , .                  , "SPy_0008 from AE004092"
    SPy_0009 , CYTOPLASM    , .                  , "SPy_0009 from AE004092"
    SPy_0010 , PSE          , tmhmm(1)           , "SPy_0010 from AE004092"
    SPy_0012 , PSE          , hmm(GW1);signalp   , "SPy_0012 from AE004092"

The first column is the SeqID which is the first token in the identifier line of the sequence in the FASTA file

The second column is the prediction, it is CYTOPLASM, MEMBRANE, PSE, or SECRETED.

The third line is a summary of features picked up by the other program:

- tmhmm(2) means 2 transmembrane helices were found
- hmm(GW1) means the GW1 motif was found in HMMER
- signalp means a secretion signal was found
- lipop means a Sp II secretion signal found with an appropriate CYS residue at the cleavage site, which will be attached to a phospholipid in the membrane

The rest of the line gives the full name as given in the rest of the identifier line in the FASTA file. This is enclosed by quotation marks "" for CSV compatability.


## Configuration file

The parameters and locations for the dependencies of _inmembrane_ are held in a file called `inmembrane.config`, which is located in the same directory as the main script. The `inmemmbrane.py` script will always look for `inmembrane.config` there.

If `inmembrane.config` is not found, running `inmembrane.py` with no arguments will generate a default `inmembrane.config`. Edit the locations of the binaries in the config file to match your system. The options are:

- the full path location of the binaries for SignalP, LipoP, TMHMM, HMMSEARCH, and MEMSAT. 
- 'organsim' indicates the type of organism - Gram+, Gram- are supported so far
- 'hmm_profiles_dir': the location of the HMM profiles for peptidoglycan binding motifs for Gram+ 
-  for HMM, you can set the cutoffs for significance, the E-value 'hmm_evalue_max', and the score 'hmm_score_min'
- the shortest length of a loop that sticks out of the peptidoglycan layer of a Gram+ bacteria. The SurfG+ determined this to be 50 amino acids for terminal loops, and twice that for internal loops, 100
- 'helix_programs' you can choose which of the transmembrane-helix prediction programs you want to use

## Installing dependencies

As it is the nature of bioinformatic programs that they are changed and updated severely with different versions, and that stable APIs and output formats are the exception rather than the norm. It is very important that you have the exact version that we have programmed against.

Required dependencies, and their versions:

- TMHMM 2.0 _or_ MEMSAT3

- SignalP 4.0

- LipoP 1.0  

- HMMER 3.0

These instructions have been tailored for Debian-based systems, in particular Ubuntu 11.10. Each of these dependencies are licensed free to academic users.

### TMHMM 2.0
Only one of TMHMM or MEMSAT3 are required, but users that want to compare transmembrane segment predictions can install both.

- Download and install TMHMM 2.0 from <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm>.

### SignalP 4.0
- Download SignalP 4.0 <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp>. You will need to fill out the form with an institutional email address and accept the academic license. The software will be emailed to you.

- Follow the installation instructions at <http://www.cbs.dtu.dk/services/doc/signalp-4.0.readme>.

### HMMER 3.0
- Download HMMER 3.0 from <http://hmmer.janelia.org/software>.
- The HMMER user guide describes how to install it. For the pre-compiled packages, this is as simple as putting the binaries on your PATH.

### LipoP 1.0
- Download LipoP 1.0 from <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?lipop>. The installation proceedure is similar to that for SignalP.

### MEMSAT3

- Download MEMSAT3 from <http://bioinfadmin.cs.ucl.ac.uk/downloads/memsat/memsat3/> (only memsat3_academic.tar is required). 

- MEMSAT3 requires NCBI BLAST ("legacy" BLAST, not BLAST+) using the SwissProt (swissprot) database.
 - Legacy BLAST can be downloaded at <ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/> installed using the instructions provided by NCBI <http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/unix_setup.html>. We have tested with version 2.2.25.
 - You will need both the 'nr' database and the 'swissprot' database, since 'swissprot' is indexed against 'nr'. (The other option is to download the FASTA version of Uniprot/Swiss-Prot from <ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz> and create your own BLAST formatted database with using the BLAST formatdb tool).

- Edit the _runmemsat_ script included with MEMSAT3 to point to the correct locations using absolute paths:
 - 'dbname' is the location of your BLAST formatted swissprot database
 - 'ncbidir' is the base directory of your BLAST installation
 - 'execdir' is the path where the MEMSAT3 executable resides
 - 'datadir' is the the path to the MEMSAT3 data directory )

Note the the 'runmemsat' script refers to PSIPRED v2, but it means MEMSAT3 - PSIPRED is not required.

## Tests

We provide a number of tests for _inmembrane_, but given the variety of different versions of the programs that _inmembrane_ depends, we are unable to comprehensively test all of them. In the _inmembrane_ directory:

     python runtest.py

