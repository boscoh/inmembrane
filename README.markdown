# inmembrane

_inmembrane_ is a pipeline for proteome annotation to predict if a protein is exposed on the surface of a bacteria. 

## Installation and Configuration

Dowload the latest version of _inmembrane_ from the github repository: <https://github.com/boscoh/inmembrane/zipball/master>. 

This package includes tests, examples, data files, docs and a few included libraries ([Beautiful Soup](http://www.crummy.com/software/BeautifulSoup/), [mechanize](http://wwwsearch.sourceforge.net/mechanize) and [twill](http://twill.idyll.org/)).

The editable parameters of _inmembrane_ are found in `inmembrane.config`, which is always located in the same directory as the main script. If no such file exists, a default `inmembrane.config` will be generated. The parameters are: 

- the path location of the binaries for SignalP, LipoP, TMHMM, HMMSEARCH, and MEMSAT. This can be the full path, or just the binary name if it is on the system path environment. Use `which` to check. 
- 'protocol' to indicate which analysis you want to use. Currently, we support:
    - `gram_pos` the analysis of surface-exposed proteins of Gram+ bacteria; 
    - `gram_neg` the sub-cellular location of Gram- bacteria.
- 'hmm_profiles_dir': the location of the HMMER profiles for any HMM peptide sequence motifs
-  for HMMER, you can set the cutoffs for significance, the E-value 'hmm_evalue_max', and the score 'hmm_score_min'
- the shortest length of a loop that sticks out of the peptidoglycan layer of a Gram+ bacteria. The SurfG+ determined this to be 50 amino acids for terminal loops, and twice that for internal loops, 100
- 'helix_programs' you can choose which of the transmembrane-helix prediction programs you want to use

We provide a number of unit tests for _inmembrane_:

    python runtest.py

As _inmembrane_ has a lot of dependencies, these tests are really useful in working out if the dependencies are installed in a way that is compatible with _inmembrane_. Since not all the binaries are needed, not all tests (and corresponding dependencies) are required for _inmembrane_ to work.

## Execution

_inmembrane_ was written in Python 2.7. It takes a FASTA input file and runs a number of external bioinformatic programs on the sequences. It then collects the output to make the final analysis, which is printed out and stored in a CSV file.  

_inmembrane_ can be run in two modes. It can be run as a command-line program:  
     
    python inmembrane.py your_fasta_file

The other way of running imembrane.py is with a custom script, such as `run_example.py` where all pertinent input is in the script itself. You can either run this on the command-line like this:

    python run_example.py

or simply double-click `run_example.py` in a file-manager. 

If you open it with a text editor, you will see that the input file and the csv file are filled in. You can change this by simply duplicating `run_example.py`, and editing the parameters in a text editor.

## Output format

The output of _inembrane_ consists of four columns of output. This is printed out on standard out and in a CSV file, which can be opened in EXCEL. If not specified, the CSV will be given the same basename as the FASTA file. The standard text output can be parsed using space delimiters (empty fields in the third column are indicated with a "."). Here's an example:

    SPy_0008   CYTOPLASM     .                   "SPy_0008 from AE004092"
    SPy_0009   CYTOPLASM     .                   "SPy_0009 from AE004092"
    SPy_0010   PSE           tmhmm(1)            "SPy_0010 from AE004092"
    SPy_0012   PSE           hmm(GW1);signalp    "SPy_0012 from AE004092"

- the first column is the SeqID which is the first token in the identifier line of the sequence in the FASTA file

- the second column is the prediction, it is CYTOPLASM, MEMBRANE, PSE, or SECRETED. A value of PSE (Potentially Surface Exposed) means that _inmembrane_ has predicted the protein to be surface exposed and will be found in a membrane-shaving experiment.

- the third line is a summary of features picked up by the other program:
    - tmhmm(2) means 2 transmembrane helices were found
    - hmm(GW1) means the GW1 motif was found in HMMER
    - signalp means a secretion signal was found
    - lipop means a Sp II secretion signal found with an appropriate CYS residue at the cleavage site, which will be attached to a phospholipid in the membrane

- the rest of the line gives the full identifier of the sequence in the FASTA file. This is enclosed by quotation marks "" in the CSV file.


## Installing dependencies

As it is the nature of bioinformatic programs that they are changed and updated severely with different versions, and that stable APIs and output formats are the exception rather than the norm. It is very important that you have the exact version that we have programmed against.

Required dependencies, and their versions:

- TMHMM 2.0 _or_ MEMSAT3
- SignalP 4.0
- LipoP 1.0 
- HMMER 3.0

These instructions have been tailored for Debian-based systems, in particular Ubuntu 11.10. Each of these dependencies are licensed free to academic users.

#### TMHMM 2.0
Only one of TMHMM or MEMSAT3 are required, but users that want to compare transmembrane segment predictions can install both.

- Download and install TMHMM 2.0 from <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm>.

#### SignalP 4.0
- Download SignalP 4.0 <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp>. You will need to fill out the form with an institutional email address and accept the academic license. The software will be emailed to you.
- Follow the installation instructions at <http://www.cbs.dtu.dk/services/doc/signalp-4.0.readme>.

#### HMMER 3.0
- Download HMMER 3.0 from <http://hmmer.janelia.org/software>.
- The HMMER user guide describes how to install it. For the pre-compiled packages, this is as simple as putting the binaries on your PATH.

#### LipoP 1.0
- Download LipoP 1.0 from <http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?lipop>. The installation proceedure is similar to that for SignalP.

#### MEMSAT3
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

## Modification guide

It is a fact of life for bioinformatics that new versions of basic tools changes output formats and API. We believe that it is an essential skill to rewrite parsers to handle the subtle but significant changes in different versions. We have written  _inmembrane_ to be easily modifiable and extensible. _Protocols_ which embody a particular high level workflow are found in `inmembrane/protocols`. 

All interaction with a specific external program or web-site have been wrapped into a single python _plugin_ module, and placed in the `plugins` directory. This contains the code to both run the program and to parse the output. We have tried to make the parsing code as concise as possible. Specifically, by using the native Python dictionary, which allows an enormous amout of flexibility, we can extract the analysis with very little code.

### __inmembrane__ development style guide:

Here are some guidelines for understanding and extending the code.

* _Confidence:_ Plugins that wrap an external program should have at least one high level test which is executed by run_tests.py. This allows new users to immediately determine if their dependencies are operating as expected.
* _Interface:_ A plugin that wraps an external program must receive a _params_ data structure (derived from `inmembrane.config`) and a _proteins_ data structure (which is a dictionary keyed by sequence id). Plugins should return a 'proteins' object.
* _Flexibility:_ Plugins should have a 'force' boolean argument that will force the analysis to re-run and overwrite output files.
* _Efficiency:_ All plugins should write an output file which is read upon invocation to avoid the analysis being re-run.
* _Documentation:_ A plugin must have a Python docstring describing what it does, what parameters it requires in the `params` dictionary and what it adds to the `proteins` data structure. See the code for examples.
* _Anal:_ Unique sequence ID strings (eg `gi|1234567`) are called 'seqid'. 'name' is ambiguous. 'prot_id' is reasonable, however conceptually a 'protein' is not the same thing as a string that represents it's 'sequence' - hence the preference for 'seqid'.
* _Anal:_ All file handles should be closed when they are no longer needed.
