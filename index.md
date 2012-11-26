## Getting it to Work

#### What is _inmembrane_?

_inmembrane_ is a program to predict if a protein is likely to be surface-exposed in a bacterial cell. It is designed to be relatively easy to install, with several different options to help you run all external binaries on your system.

#### How can we help you?

We envisage _inmembrane_ being of use in three ways:

1. You want to do a quirk and dirty annotation of your proteome for cell-surface exposure in order to cross-check your experiment. Keep reading on this page.

2. You want to play around with different ways of cell localization assignment. Check out our parameterization guide.

3. You want to make your own bioinformatic work-flow, and would like to use _inmembrane_ as a useful Python starting point. Have a look at our programming model, paper, and API guide.


#### Combining analyses into a workflow

The cell-surface exposure algorithm is taken from the paper (Barinove et al. 2009. Proteomics 9:61-73). It is a straightforward case of a simple work-flow defined over a number of external binaries that do:

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

It should look something like this: 

    SPy_0008  CYTOPLASM(non-PSE)  .                         SPy_0008 from AE004092
    SPy_0010  PSE-Membrane        tmhmm(1)                  SPy_0010 from AE004092
    SPy_0012  PSE-Cellwall        hmm(GW2|GW3|GW1);signalp  SPy_0012 from AE004092
    SPy_0016  MEMBRANE(non-PSE)   tmhmm(12)                 SPy_0016 from AE004092
    SPy_0019  SECRETED            signalp                   SPy_0019 from AE004092

There are 4 columns to the result: 

1. The sequence identifier (seqid) is the first word in the id line of the FASTA file
2. The cell-localization category:
    - PSE: means Potentially Surface Exposed, this can arise from
        - transmembrane helical protein with long longs
        - lipoprotein attached to the lipid layer
        - protein with motif that attaches to the glycan cell-wall
    - SECRETED: contains a secretion signal
    - CYTOPLASM: basically none of the above
3. Summary of the results from the external bioinformatic programs. In this case, we have hits from `tmhmm`, `hmm`, and `signalp`
4. The description of the protein taken from the FASTA file.

