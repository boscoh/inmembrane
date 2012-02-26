# _inmembrane_

_inmembrane_ is a pipeline for proteome annotation to predict if a protein  is exposed on the surface of a bacteria. It takes a FASTA input file and runs  a number of external bioinformatic analysis programs on the sequences. It then collects the output to make the final analysis.

The program is written to use as much of the standard Python library as possible, and specifically uses native data structures which provide an incredible amount of flexibility.

In particular, the goal of _inmembrane_ is to write the parsing code for each program as concise as possible. It is better to write a new parser for different versions of a program than to to try to write a monolithic parser that handles many versions. 

## Installation

_inmembrane_ is written in Python, and runs in Python 2.7. It can be run in two modes:

1. as a command-line program in the _inmembrane_ directory: python _inmembrane_.py 
2. through a wrapped python program, that can be double-clicked in a file-manager on many systems

_inmembrane_ is basically a wrapper around a large number of standard bioinformatics tools. 

The parameters and locations for these dependencies are held in a file called _inmembrane_.config, which is located in the same directory as the main script. 
If _inmembrane_.config is not found, running _inmembrane_ with no arguments will generate a default _inmembrane_.config file. Edit the locations of the binaries in the config file to match your system.

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

## Tests

We provide a number of tests for _inmembrane_, but given the variety of different versions of the programs that _inmembrane_ depends, we are unable to comprehensively test all of them. In the _inmembrane_ directory:

     python runtest.py

