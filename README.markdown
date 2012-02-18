# INMEMBRANE

Inmembrane.py is a program that is designed to bind the results of various other bioinformatic programs together to form an analysis.

The philosophy of inmembrane.py is that the standard Python library provides a clean class of functions to process the relatively straightforward analysis and that simple dictionaries and python functions should be able to provide all the analysis required for the analysis.

the current set of scripting languages provide such a clean interface to the tasks required for scripting that we believe that it represents a significant obscuratuon to wrap these commonly used fuctions in further abstraction. 

You cannot any better expressability than in the native language constructs. 

As well experience has shown that programs can and most obviously will change over iterative versions. it is better to accept that output will change and allow easy additions to the parser rather than try to write a monolithic parser that will premot all future changes. light and adaptability. 


## Installation

inmembrane.py is a plain Python script, and has been tested with Python 2.7. 

Since inmembrane is basically a wrapper, it requires a lot programs to be installed. The parameters and locations for these dependencies are held in a file called inmembrane.config, which is located in the same directory as the main script. If this is not found, simply running INMEMBRANE with no arguments will generate the default inmembrane.config file. Edit the file to match your system after this file is generated.

unit tests to see if your bioinformatics works with inmembrane. 

run a test with a short example and see if it crashes. 

if it cannot find the binary, it will quit and inform you of this. 

One problem with testing something like imembrane is that there is a wide variety of ways in which its dependencies can be installed. For instance, memsat3 in turn runs a local copy of NCBI Blast , it is difficult to test because memsat3 uses blast which in turn depends on the exact version of the blast binaries and databases that you install. 

## Dependencines

As it is the nature of bioinformatic programs that they are changed and updated severely with different versions, and that stable APIs and output formats are the exception rather than the norm. It is very important that you have the exact version that we have programmed against.

Nevertheless, we believe that it is essential skill to rewrite parsers to handle the subtle but significant changes in different versions. 

Required dependencies, and their versions:

- TMHMM 2.0 _or_ MEMSAT3

- SignalP 4.0

- LipoP 1.0  

- HMMER 3.0

## Installing dependencies

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
