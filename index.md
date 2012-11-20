
#### What is _inmembrane_?

_inmembrane_ is a program to predict if a protein is likely to be surface-exposed in a bacterial cell. It is designed to be relatively easy to install, with several different options to help you run all external binaries on your system.

#### How can we help you?

We envisage _inmembrane_ being of use in three ways:

1. You want to do a quirk and dirty annotation of your proteome for cell-surface exposure in order to cross-check your experiment. Keep reading on this page.

2. You want to play around with different ways of cell localization assignment. Check out our parameterization guide.

3. You want to make your own bioinformatic work-flow, and would like to use _inmembrane_ as a useful Python starting point. Have a look at our programming model, paper, and API guide.


#### Combining analyses into a workflow

The cell-surface exposure algorithm is taken from the paper (). It is a straightforward case of a simple work-flow defined over a number of external binaries that do:

- transmembrane helix prediction,
- lipoprotein prediction,
- secretion signal prediction
- cell-wall binding motif prediction

Although there exists several generic bioinformatic workflow packages, the deduction of cell-surface exposure from the results of these programs is a little more involved. Here, we have written this in a Python package, and in a manner that is extensible, and potentially repackageable for other bioinformatic analysis. 

#### Ease of use

We have specifically chosen Python for the overall ease of use, and readability. Read our paper for more details. One of the most difficult aspects of installing a workflow like this is the installation of multiple external binaries. We have made a large effort to incorporate web-based modules, which in most cases, should significantly simplify the process of installation.

#### How to install?

Well, to install the program, you'll have to have install

 - Python: it should be on most systems, but if not, go to <http://www.python.org>
 - pip: <http://pypi.python.org/pypi/pip>
 - HMMER: this is the only bioinformatic software that definitely needs to be installed <http://hmmer.janelia.org/software>

Once these are installed on your system, you should be able to do this to install _inmembrane_:

    >>> pip install inmembrane

#### Running _inmembrane_

Let's say you have all the protein sequences you want to analyze in a FASTA text file called `proteome.fasta`.

Inmembrane can be run from the command line, assuming you have an internet connection, and all the web-sites _inmembrane_ relies on are up:

    >>> inmembrane_scan proteome.fasta

The results will be written to proteome.csv, which can be opened by Microsoft Excel.

#### Interpreting output

The program is to predict the 

#### The parameters in the program

Here are the parameters in the program, as indicated in the `inmembrane.config` file and in the get_params() function of the API. 

_general parameters_

- 'fasta': the pathname of the fasta file holding the sequences
- 'csv': the name of the output file in CSV (Excel-compatible) format. Defaults to a file using the same filename as 'fasta'
- 'out_dir': the directory where all intermediate output files are located
- 'protocol': switches between different protocols for analysising the proteomes

_optional binaries_

- 'signalp4_bin', 'lipop1_bin', 'tmhmm_bin', 'memsat3_bin': the full pathname of the binaries used for the analysis. If empty, it is not run, and the web version is used instead

_transmembrane helix predictors_

- 'helix_programs': a choice of which transmembrane-helix predictors are used. Currently, the program loops over all designated helix predictors and chooses the longest predicted loops.

_key parameters to determine surface-exposure_

- 'terminal_exposed_loop_min': this is the length in terms of amino acids given to define an external loop of a protein to be long enough to stick out of the bacterial cell-wall to be considered surface-epxosed
- 'internal_exposed_loop_min': similarly to above, except is double the length as internal loops need to go up, and down again, for continuity with the rest of the protein.

_HMMER parameters for signal proteins for cell-wall attachment_

- 'hmmsearch3_bin': the binary for HMMER
- 'hmm_profiles_dir': the directory that holds the HMMER profiles for cell-wall signals
- 'hmm_evalue_max': the maximum expectation value cutoff for acceptable prediction in HMMER
- 'hmm_score_min': the minimum score cutoff for acceptatble prediction in HMMER

&beta;_-barrel predictors_

- 'barrel_programs': list of programs to do &beta-barel ['bomp', 'tmbetadisc-rbf'],
- 'bomp_clearly_cutoff': 3, # >= to this, always classify as an OM(barrel)
- 'bomp_maybe_cutoff': 1, # must also have a signal peptide to be OM(barrel)
- 'tmbhunt_clearly_cutoff': 0.95,
- 'tmbhunt_maybe_cutoff': 0.5,
- 'tmbetadisc_rbf_method': 'aadp', # aa, dp, aadp or pssm

#### Advanced Users

