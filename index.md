Need a change.

#### Problem

You have a bunch of protein sequences in fasta form. Your organism has not been well studied, and you'd like to see if any of them are likely to be exposed on the cell surface. 

#### Combining analyses into a workflow

Actually, it's quite easy to predict surface-exposure by combining the results of several _other_ prediction programs. However, to do this will require some kind of work-flow to combine all the results of the required prediction programs, such as helix predictors and lipo-protein predictiors.

There are several general approaches to develop bioinformatic workflows for protein sequence analysis, many of them requiring a high-level of expertise for installation. Here, we provide a Python-based solution that relies on web-based services, that we believe represents an easy, transparent, and extensible approach.

#### Ease of use

We've taken the pain of installing all the external binaries by hooking into the web-based versions of these programs. Installing local binaries is painful, and sensitive to platform inconsistencies. By focusing on web-basesd services, we can eliminate a lot of this pain.

#### Why would you use this?

We see primarilty three types of users:

1) You want to do a quirk and dirty annotation of your proteome and get a prediction of cell-surface exposure. Perhaps this is a run-once only script.

2) You want to play around with different ways of cell localization assignment. 

3) You want to make your own bioinformatic work-flow, and would like a useful starting point.

#### How to install?

Well, to install the program, you'll have to have 

* Python
* pip: 
* HMMER: make sure this is available on the command-line

Once these are installed on your system, to install _inmembrane_, you should be able to do this:

    >>> pip install inmembrane

Once installed, inmembrane can be run from the command line, assuming you have an internet connection:

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

