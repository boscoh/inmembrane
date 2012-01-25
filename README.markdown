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

The programs that we requre, and their versions:

- MEMSAT3 <http://bioinf.cs.ucl.ac.uk/?id=756> Edit the runmemsat script included with MEMSAT3 to point to the correct locations with absolute paths (eg dbname, ncbidir, execdir, datadir)
  - MEMSAT3 in turns requires the NCBI BLAST using the SwissProt database

- SignalP 1.0

- LipoP 1.0  

- HMMER 3.0

- TMHMM 2.0