inmembrane
==========

*inmembrane* is a pipeline for proteome annotation to predict if a
protein is exposed on the surface of a bacteria. It orchestrates the 
analysis of protein sequences to provides a summary of which targets may 
be surface exposed based on predicted subcellular localization signals and 
membrane topology. Currently protocols have been implemented for gram+ and
gram- bacterial proteomes.

Typical usage is via the script ``inmembrane_scan.py``, eg::

    $ inmembrane_scan.py mysequences.fasta


The provided sequences 
(in `FASTA format <http://en.wikipedia.org/wiki/FASTA_format>`_) 
are subjected to an number of sequence analyses using external
programs (see below) and the result summarized like::

  SPy_0008  CYTOPLASM(non-PSE)  .                         SPy_0008 from AE004092
  SPy_0010  PSE-Membrane        tmhmm(1)                  SPy_0010 from AE004092
  SPy_0012  PSE-Cellwall        hmm(GW2|GW3|GW1);signalp  SPy_0012 from AE004092
  SPy_0016  MEMBRANE(non-PSE)   tmhmm(12)                 SPy_0016 from AE004092
  SPy_0019  SECRETED            signalp                   SPy_0019 from AE004092


As well as output to stdout, this will generate a summary CSV file 
``mysequences.csv``and a directory ``mysequences`` containing output
files generated during the run.

Although *inmembrane* is primarily designed to be used as a stand alone
program, it can also be used as a library like::

  import inmembrane
  params = inmembrane.get_params()
  params['fasta'] = "input.fasta"
  annotations = inmembrane.process(params)

where ``annotations`` is a dictionary of the results, with protein sequence IDs as
keys.

You can also test the functionality of the analysis plugins
that are part of *inmembrane* by typing::

    $ inmembrane_scan.py --test

This can be useful for determining which binary dependences
are correctly installed, or exposing any broken / offline web services
required for a particular analysis.

Installation and Configuration
------------------------------

The latest stable release of *inmembrane* can be installed via 
pip, or the bleeding edge from Github.

Via pip::

    $ sudo pip install inmembrane

Or from Github::

    $ git clone http://github.com/boscoh/inmembrane.git
    $ cd inmembrane
    $ sudo python setup.py install

The package includes tests, examples, data files, docs.
HMMER3 is the only *required* external dependency, however
for large analyses (multiple proteomes) it is suggested 
that local versions of other analysis programs are installed 
rather than relying on web services (see Installing dependencies_ below).

The editable parameters of *inmembrane* are found in
``inmembrane.config``, which is always located in the same
directory as the main script. If no such file exists, a default
``inmembrane.config`` will be generated. By default, you probably
don't need to change anything.

The parameters are:

-  the path location of the binaries for SignalP, LipoP, TMHMM,
   HMMSEARCH, and MEMSAT. This can be the full path, or just the
   binary name if it is on the system path environment. Use ``which``
   to check.
-  'protocol' to indicate which analysis you want to use.
   Currently, we support:
   
   -  ``gram_pos`` the analysis of surface-exposed proteins of Gram+
      bacteria;
   -  ``gram_neg`` annotation of subcellular localization and inner
      membrane topology classification for Gram- bacteria

-  'hmm\_profiles\_dir': the location of the HMMER profiles for any
   HMM peptide sequence motifs
-  for HMMER, you can set the cutoffs for significance, the E-value
   'hmm\_evalue\_max', and the score 'hmm\_score\_min'
-  the shortest length of a loop that sticks out of the
   peptidoglycan layer of a Gram+ bacteria. The SurfG+ determined this
   to be 50 amino acids for terminal loops, and twice that for
   internal loops, 100
-  'helix\_programs' you can choose which of the
   transmembrane-helix prediction programs you want to use

Output format
-------------

The output of *inmembrane* ``gram_pos`` protocol consists of four
columns of output. This is printed to stdout and written as a CSV
file, which can be opened in spreadsheet software such as EXCEL.
The standard text output can be parsed using space delimiters
(empty fields in the third column are indicated with a ".").
Logging information are prefaced by a '#' character, and is sent to
stderr.

Here's an example::

  SPy_0008  CYTOPLASM(non-PSE)  .                         SPy_0008 from AE004092
  SPy_0009  CYTOPLASM(non-PSE)  .                         SPy_0009 from AE004092
  SPy_0010  PSE-Membrane        tmhmm(1)                  SPy_0010 from AE004092
  SPy_0012  PSE-Cellwall        hmm(GW2|GW3|GW1);signalp  SPy_0012 from AE004092
  SPy_0013  PSE-Membrane        tmhmm(1)                  SPy_0013 from AE004092
  SPy_0015  PSE-Membrane        tmhmm(2)                  SPy_0015 from AE004092
  SPy_0016  MEMBRANE(non-PSE)   tmhmm(12)                 SPy_0016 from AE004092
  SPy_0019  SECRETED            signalp                   SPy_0019 from AE004092


-  the first column is the SeqID which is the first token in the
   identifier line of the sequence in the FASTA file

-  the second column is the prediction, it is CYTOPLASM(non-PSE),
   MEMBRANE(non-PSE), PSE-Cellwall, PSE-Membrane, PSE-Lipoprotein or
   SECRETED. Any 'PSE' (Potentially Surface Exposed) annotation means
   that based on the predicted topology, the protein is likely to be
   surface exposed and will be protease accessible in a
   membrane-shaving experiment.

-  the third line is a summary of features detected by external
   tools:
   
   -  tmhmm(2) means 2 transmembrane helices were found by TMHMM
   -  hmm(GW2\|GW3\|GW1) means that the GW1, GW2 and GW3 motifs were
      found by HMMER hmmsearch
   -  signalp means a secretion signal was found SignalP
   -  lipop means a Sp II secretion signal found by LipoP with an
      appropriate CYS residue at the cleavage site, which will be
      attached to a phospholipid in the membrane

-  the rest of the line gives the full identifier of the sequence
   in the FASTA file.


Installing dependencies
-----------------------
.. _dependencies:

While *inmembrane* only requires a local installation of HMMER 3.0
and can used web services for TMHMM, SignalP, LipoP and various
OMP beta-barrel predictors, for large scale analyses (5000 sequences+)
it is suggested that locally installed versions are used in the interest
of speed, at to be polite to publically available web services.

As it is the nature of bioinformatic programs that they are changed
and updated severely with different versions, stable APIs with
consistent output formats are the exception rather than the norm.
It is very important that you have the exact version that *inmembrane*
is written to interoperate with.

Required dependencies, and their versions:


-  TMHMM 2.0 *or* MEMSAT3
-  SignalP 4.0
-  LipoP 1.0
-  HMMER 3.0

These instructions have been tailored for Debian-based systems, in
particular Ubuntu 11.10. Each of these dependencies are licensed
free to academic users.

TMHMM 2.0
^^^^^^^^^

Only one of TMHMM or MEMSAT3 are required, but users that want to
compare transmembrane segment predictions can install both.


-  Download and install TMHMM 2.0 from
   http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm.

SignalP 4.0
^^^^^^^^^^^


-  Download SignalP 4.0
   http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?signalp. You will need
   to fill out the form with an institutional email address and accept
   the academic license. The software will be emailed to you.
-  Follow the installation instructions at
   http://www.cbs.dtu.dk/services/doc/signalp-4.0.readme.

HMMER 3.0
^^^^^^^^^


-  Download HMMER 3.0 from http://hmmer.janelia.org/software.
-  The HMMER user guide describes how to install it. For the
   pre-compiled packages, this is as simple as putting the binaries on
   your PATH.

LipoP 1.0
^^^^^^^^^


-  Download LipoP 1.0 from
   http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?lipop. The
   installation proceedure is similar to that for SignalP.

MEMSAT3
^^^^^^^


-  Download MEMSAT3 from
   http://bioinfadmin.cs.ucl.ac.uk/downloads/memsat/memsat3/ (only
   memsat3\_academic.tar is required).
-  MEMSAT3 requires NCBI BLAST ("legacy" BLAST, not BLAST+) using
   the SwissProt (swissprot) database.
-  Legacy BLAST can be downloaded at
   ftp://ftp.ncbi.nlm.nih.gov/blast/executables/release/LATEST/
   installed using the instructions provided by NCBI
   http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/unix_setup.html. We
   have tested with version 2.2.25.
-  You will need both the 'nr' database and the 'swissprot'
   database, since 'swissprot' is indexed against 'nr'. (The other
   option is to download the FASTA version of Uniprot/Swiss-Prot from
   ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
   and create your own BLAST formatted database with using the BLAST
   formatdb tool).

-  Edit the *runmemsat* script included with MEMSAT3 to point to
   the correct locations using absolute paths:
-  'dbname' is the location of your BLAST formatted swissprot
   database
-  'ncbidir' is the base directory of your BLAST installation
-  'execdir' is the path where the MEMSAT3 executable resides
-  'datadir' is the the path to the MEMSAT3 data directory )


(Note the the 'runmemsat' script refers to PSIPRED v2, but it means
MEMSAT3 - PSIPRED is NOT required).

Python libraries
^^^^^^^^^^^^^^^^

*inmembrane* depends on the following Python libraries (
`Beautiful Soup <http://www.crummy.com/software/BeautifulSoup/>`_,
`mechanize <http://wwwsearch.sourceforge.net/mechanize>`_ and
`twill <http://twill.idyll.org/>`_) and
`Suds <https://fedorahosted.org/suds/>`_ ). Pip should handle
installing these for you automatically.

Modification guide
------------------

It is a fact of life for bioinformatics that new versions of basic
tools changes output formats and API. We believe that it is an
essential skill to rewrite parsers to handle the subtle but
significant changes in different versions. We have written
*inmembrane* to be easily modifiable and extensible. *Protocols*
which embody a particular high level workflow are found in
``inmembrane/protocols``.

All interaction with a specific external programs or web services have
been wrapped into a single python *plugin* module, and placed in
the ``inmembrane/plugins`` directory. This contains the code to both run the
program and to parse the output. We have tried to make the parsing
code as concise as possible. Specifically, by using the native
Python dictionary, which allows an enormous amout of flexibility,
we can collate the results of various analyses with very little code.

**inmembrane** development style guide:
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here are some guidelines for understanding and extending the code.


-  *Confidence:* Plugins that wrap an external program should have
   at least one high level test which is executed by run\_tests.py.
   This allows new users to immediately determine if their
   dependencies are operating as expected.
-  *Interface:* A plugin that wraps an external program must
   receive a *params* data structure (derived from
   ``inmembrane.config``) and a *proteins* data structure (which is a
   dictionary keyed by sequence id). Plugins should return a
   'proteins' object.
-  *Flexibility:* Plugins should have a 'force' boolean argument
   that will force the analysis to re-run and overwrite output files.
-  *Efficiency:* All plugins should write an output file which is
   read upon invocation to avoid the analysis being re-run.
-  *Documentation:* A plugin must have a Python docstring
   describing what it does, what parameters it requires in the
   ``params`` dictionary and what it adds to the ``proteins`` data
   structure. See the code for examples.
-  *Anal:* Unique sequence ID strings (eg ``gi|1234567``) are
   called 'seqid'. 'name' is ambiguous. 'prot\_id' is reasonable,
   however conceptually a 'protein' is not the same thing as a string
   that represents it's 'sequence' - hence the preference for 'seqid'.
-  *Anal:* All file handles should be closed when they are no
   longer needed.


