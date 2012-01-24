Inmembrane.py is a program that is designed to bind the results of various other bioinformatic programs
together to form an analysis.

The philosophy of inmembrane.py is that the standard Python library provides a clean class of functions
to process the relatively straightforward analysis and that simple dictionaries and python functions
should be able to provide all the analysis required for the analysis.

# Installation

inmembrane.py is a plain Python script, and has been tested with Python 2.7. 

Since inmembrane is basically a wrapper, it requires a lot programs to be installed
The parameters and locations for these dependencies are held in a file
called inmembrane.config, which is located in the same directory as the main script.

If this is not found, simply running INMEMBRANE with no arguments will generate
the default inmembrane.config file. Edit the file to match your system after this
file is generated.

# Dependencines

As it is the nature of bioinformatic programs that they are changed and updated
severely with different versions, and that stable APIs and output formats are
the exception rather than the norm. It is very important that you have the exact
version that we have programmed against.

Nevertheless, we believe that it is essential skill to rewrite parsers to handle
the subtle but significant changes in different versions. 

The programs that we requre, and their versions:

* MEMSAT3 ( http://bioinf.cs.ucl.ac.uk/?id=756 
            http://bioinfadmin.cs.ucl.ac.uk/downloads/memsat/memsat3/)
  Edit the runmemsat script included with MEMSAT3
  to point to the correct locations with absolute paths
  (eg dbname, ncbidir, execdir, datadir)

  - MEMSAT3 in turns requires the NCBI BLAST using the SwissProt database

* SignalP 1.0

* LipoP 1.0 - 

* HMMER 3.0

* TMHMM 2.0