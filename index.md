## What is _inmembrane_?

_inmembrane_ is a program to predict if a set of proteins are likely to be surface-exposed in a bacterial cell based on their predicted subcellular localization and membrane topology. It is designed to be relatively easy to install, with options for running the analysis using web services and/or locally installed sequence analysis packages.

_inmembrane_ runs on Linux, OS X and other Unix-based systems. It has been tested on Ubuntu Linux (11.04+) and Mac OS X (10.6+). It hasn't been tested on Windows, and is unlikely to run without changes.

### Quickstart for the bold & impatient

    $ sudo pip install inmembrane
    $ inmembrane_scan --test
    $ inmembrane_scan my_gram_pos_proteome.fasta

<em>(Some tests may fail if you don't have various external binaries installed, but don't despair - the only external program **required** by is HMMER v3, everything else uses a web service in the default configuration)</em>

### How can inmembrane help you?

We envisage _inmembrane_ being of use in three ways:

1. Quirk and dirty annotation of a proteome in order to cross-check your experiment (eg, gaining an overview of hits from a gram-positive bacterial cell-surface shaving experiment).

2. You want to tweak your own rule-based method for subcellular localization assignment from sequence. Check out our <a href="customization.html">parameterization guide</a>.

3. You want to make your own protein sequence analysis work-flow, and would like to use _inmembrane_ as a useful Python starting point. Have a look at our [programming model and API guide](api.html).


### Combining analyses into a workflow

The cell-surface exposure algorithm used by the _inmembrane_ 'gram-positive' protocol is based on the method described in <a href="http://dx.doi.org/10.1002/pmic.200800195">(Barinov et al. 2009. Proteomics 9:61-73)</a>. The workflow involves annotating the provided protein sequences using a number of external programs and/or web services providing:

- transmembrane helix prediction (eg TMHMM)
- lipoprotein prediction (eg LipoP)
- secretion signal prediction (eg SignalP)
- cell-wall binding motif prediction (eg via HMMER)

_inmembrane_ simplifies collating the output and presenting an overview which allows biologists to immediately gain insight into the likely localization and surface exposure of their bacterial proteome of interest.

## Using _inmembrane_

### Installation

Installing and running _inmembrane_ requires:

 - Python: it should be on most modern systems, but if not, go to <http://www.python.org>
 - pip: <http://pypi.python.org/pypi/pip>
 - HMMER: this is the only sequence analysis software that is required to be locally installed <http://hmmer.janelia.org/software>

On Debian/Ubuntu systems, you can install the minimal dependencies using:

    $ sudo apt-get install python-pip hmmer

Once these are installed on your system, _inmembrane_ can be installed with:

    $ sudo pip install inmembrane

### Running _inmembrane_

Let's say you have all the protein sequences you want to analyze in a FASTA text file called `my_proteome.fasta`.

_Inmembrane_ is run from the command line like:

    $ inmembrane_scan my_proteome.fasta

Assuming the web services _inmembrane_ relies on are available, the analysis should complete within seconds or minutes.

A summary is printed to the terminal, and the results will also be written to `my_proteome.csv`, which can be opened by a spreadsheet application such as Microsoft Excel. Citations for all external programs and services are written to `my_proteome/citations.txt` - this is very handy and important when it comes to writing up a publication using the results of _inmembrane_ and the software it orchestrates.

**Pro-tip:** The directory created by _inmembrane_ (eg named `my_proteome`), as well as containing the `citations.txt` list, also contains cached output files from the individual external analyses, the original configuration file and the input FASTA file. If you run _inmembrane_ a second time on the same FASTA file, the results will be read from this directory rather than re-running each external analysis program. To re-run the analysis from scratch (with say, different parameters), you must remove (or rename) this directory. To re-run just part of the analysis you can remove individual files from this directory.

### Interpreting the output

It should look something like this: 

    SPy_0008  CYTOPLASM(non-PSE)  .                         SPy_0008 from AE004092
    SPy_0010  PSE-Membrane        tmhmm(1)                  SPy_0010 from AE004092
    SPy_0012  PSE-Cellwall        hmm(GW2|GW3|GW1);signalp  SPy_0012 from AE004092
    SPy_0016  MEMBRANE(non-PSE)   tmhmm(12)                 SPy_0016 from AE004092
    SPy_0019  SECRETED            signalp                   SPy_0019 from AE004092

There are 4 columns to the result: 

1. The sequence identifier (SeqID) is the first word in the id line of the FASTA file
2. The cell-localization category:
    - PSE: means Potentially Surface Exposed, this can arise from
        - transmembrane helical protein with long extracellular loops
        - lipoprotein attached to the membrane
        - protein with motif that associates with the glycan cell-wall
    - SECRETED: contains a secretion signal, but not the features above
    - CYTOPLASM: contains none of the above features
3. Summary of the results from the external bioinformatic programs. In this case, we have hits from `tmhmm`, `hmm`, and `signalp`
4. The description of the protein taken from the FASTA file.

