# inmembrane

_inmembrane_ is a pipeline for proteome annotation to predict if a protein is exposed on the surface of a bacteria. 

As _inmembrane_ is a glue script that runs a lot of external dependencies, we have tried to make the parsing code for each program as concise as possible. The program is written to use as much of the standard Python library as possible, and specifically uses native data structures which provide an incredible amount of flexibility.

It is a fact of life for bioinformatics that new versions of basic tools changes output formats and API. We believe that it is an essential skill to rewrite parsers to handle the subtle but significant changes in different versions. Hopefully _inmembrane_ makes this as simple as possible.


## Execution

_inmembrane_ was written in Python 2.7. It takes a FASTA input file and runs  a number of external bioinformatic analysis programs on the sequences. It then collects the output to make the final analysis. It can be run in two modes. 

It can be run as a command-line program in the _inmembrane_ directory:  
     
    python inmembrane.py your_fasta_file

It can also be run as an auxiliary script that can be double-clicked in a file-manager. An example is give in:

    run_example.py


## Configuration file

The parameters and locations for the dependencies of _inmembrane_ are held in a file called `inmembrane.config`, which is located in the same directory as the main script. The `inmemmbrane.py` script will always look for `inmembrane.config` there.

If `inmembrane.config` is not found, running `inmembrane.py` with no arguments will generate a default `inmembrane.config`. Edit the locations of the binaries in the config file to match your system.


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

