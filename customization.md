# Configuration of _inmembrane_ 

By default, _inmembrane_ will run the Gram+ protocol using remote web services for the analysis (and HMMER locally) without requiring any additional configuration.

On the first use of the program, a default `inmembrane.config` configuration file will be generated. Once created, you can edit this `inmembrane.config`. Subsequent use of the program will reference the existing `inmembrane.config`. This file contains all the user adjustable parameters (as a Python dictionary). To run the Gram- protocol, use locally installed copies of SignalP/LipoP/TMHMM or tweak the external loop threshold, you may want to edit the configuration file to your needs.

The default parameters are as following, and they will be explained in sections below:

    {
      'fasta': '',
      'csv': '',
      'out_dir': '',
      'protocol': 'gram_pos', # 'gram_neg'
      
    #### Signal peptide and transmembrane helix prediction
      'signalp4_bin': '',
    # 'signalp4_bin': 'signalp',
      'lipop1_bin': '',
    # 'lipop1_bin': 'LipoP',
      'tmhmm_bin': '',
    #  'tmhmm_bin': 'tmhmm',
      'memsat3_bin': 'runmemsat',
      'helix_programs': ['tmhmm'],
    # 'helix_programs': ['tmhmm', 'memsat3'],
      'terminal_exposed_loop_min': 50, # unused in gram_neg protocol
      'internal_exposed_loop_min': 100, # try 30 for gram_neg
      
    #### Sequence similarity and motif prediction
      'hmmsearch3_bin': 'hmmsearch',
      'hmm_profiles_dir': '%(hmm_profiles)s',
      'hmm_evalue_max': 0.1,
      'hmm_score_min': 10,
      
    #### Outer membrane beta-barrel predictors
      'barrel_programs': ['bomp', 'tmbetadisc-rbf'],
    # 'barrel_programs': ['bomp', 'tmbetadisc-rbf', 'tmbeta'],
      'bomp_clearly_cutoff': 3, # >= to this, always classify as an OM(barrel)
      'bomp_maybe_cutoff': 1, # must also have a signal peptide to be OM(barrel)
      'tmbetadisc_rbf_method': 'aadp', # aa, dp, aadp or pssm
    }

Existing `inmembrane.config` files can be changed. Since this file is in fact just Python code, normal Python syntax rules apply. Everything after a hash (#) character on a line is ignored - these are used as comments, and to used indicate alternative configurations that you might want to use.

The parameters will be explained below, grouped in terms of related parameters.

### General Parameters

The input to _inmembrane_ is a FASTA formatted set of protein sequences. There are currently two different protocols (or 'workflows') - one for Gram+ bacteria (the default), and one for Gram- bacteria. The results are saved to a CSV file as an output.

- `fasta`: the pathname of the FASTA file holding the protein sequences
- `protocol`: switches between different protocols for analysing the proteomes, currently you can choose between `gram_pos` and `gram_neg`.
- `csv`: the name of the output file in CSV (Excel-compatible) format. Defaults to a file using the same base filename as 'fasta'

Both `fasta` and `csv` should usually be left blank by default, since the input file is usually specified directly on the commandline.

For potential debugging, it's important that all intermediate files are saved, and they are saved in: 

- `out_dir`: the directory where all intermediate output files are located. Defaults to the same filename as `fasta` (without the .fasta extension)

### Optional - Locally installed binaries or remote web services

One of the problems that _inmembrane_ tackles is the complexity of using and installing multiple pieces of sequence analysis software. For older programs, the only option was to install the local binaires of all the external programs. This raises problems involved with getting the binary (or source) for the right platform, compilation with the right tools, and installation and configuration in the correct location. In contrast, many of these tools now provide a web-interface. In _inmembrane_, we have added facilities to use remote web services, when provided, whether through a RESTful interface, or through a SOAP interface.

By default, _inmembrane_ will use web services, specified by setting the `*_bin` parameter as empty string. To use a local binary binary version of a program, the appropriate `*_bin` parameter must be set to the name (or full path) of the executable.

- `signalp4_bin`, `lipop1_bin`, `tmhmm_bin`, `memsat3_bin`: the full pathname of the binaries used for the analysis. If empty, the web version is used instead.

### The Potentially Surface Exposed (PSE) Algorithm

The key part of the surface-exposure prediction evaluates the length of loops in the extramembrane parts of membrane proteins. This analysis makes most sense for Gram+ bacteria with a single membrane, where any membrane protein can potentially expose loops to the extracellular environment. This evaluation requires the key parameter of loop cutoff that determines the minimum length of a loop that is required to potentially pass through the glycan layer. 

In the original _SurfG+_ paper [(Barinov et al. 2009. Proteomics 9:61-73)](http://dx.doi.org/10.1002/pmic.200800195), it was found that a linear polypeptide length of ~50 amino acids most closely matched the length able to traverse the cell wall and become succeptible to protease shaving in the model Gram+ bacterium studied. Therefore for internal loops between transmembrane segments, this gives an accessible length threshold of ~100 amino acids, which allows the loop to pass through the cell-wall layer and back down to the membrane again.

- `terminal_exposed_loop_min`: this is the length (in amino acids) given to define an N- or C- terminal segment of a protein that is long enough to stick out of the bacterial cell-wall to be considered surface-exposed
- `internal_exposed_loop_min`: the length of a loop between two membrane embedded segments able to traverse the cell wall and present as surface exposed

These values should be modified based on the bacterial species (based on cell wall thickness) and type of experiment being used to determine surface exposure (protease shaving vs. epitope mapping vs. chemical crosslinking).

There are several different ways to identify proteins with a subcellular localization that may display loops that protruding outside the glycan cell-wall of the bacteria. These will be discussed in the following sections.

### Proteins attached to the Glycan Cell-wall and the Lipid bilayer

The clearest example of proteins that may protrude outside of the Gram+ glycan layer are proteins that possess known cell-wall attachment motifs. Motifs for domains known to bind to the glycan cell-wall are taken from PFAM as HMMER profiles, and are included in the _inmembrane_ installation. For the Gram+ protocol, these are: CW_binding_1, Gram_pos_anchor, GW2, LPxTG, LysM, PG_binding_1, PG_binding_3, CW_binding_2, GW1, GW3, LPxTG_PS50847, NLPC_P60, PG_binding_2 and SLH.

HMMER is used to search the profiles against the input query sequences. The cutoffs for the identification of a motif are given by `hmm_evalue_max` and `hmm_score_min`.

- `hmmsearch3_bin`: the binary for HMMER
- `hmm_evalue_max`: the maximum expectation value cutoff for acceptable prediction in HMMER
- `hmm_score_min`: the minimum score cutoff for acceptatble prediction in HMMER

Another class of proteins associated with the membrane surface membrane are lipoproteins, which are posttranlationally modified to include a covalently attached N-terminal lipid moeity. LipoP is used to identify the N-terminal lipidation motif, found after a secretion signal. 

- `lipop_bin`: the location of the binary of LipoP. If this is given as '' (default), then _inmembrane_ uses the web service at <http://www.cbs.dtu.dk/services/LipoP/>.

### Transmembrane-&alpha;-helix Predictors

The other main class of protein in Gram+ proteins with extracellular loops are integral membrane proteins containing transmembrane &alpha;-helices. The loops between transmembrane segments, as well as N-terminal and C-terminal non-membrane embedded regions are all candidates to be exposed outside the glycan layer.

There is no _de facto_ standard transmembrane helix predictor, and thus _inmembrane_ has been written to allow pluggable options for transmembrane helix prediction. To allow for this, _inmembrane_ allows the results of more than one helix predictor, and takes a greedy approach to annotating potentially surface exposed proteins - the prediction that gives the longest extracellular loops will be used. The parameter 'helix_programs' lists all the transmembrane helix predictors to be used.

- `helix_programs`: a list of transmembrane-helix predictors to be used. Currently, _inmembrane_ can handle `tmhmm` and `memsat3`.

### &beta;-barrel predictors
These parameters only apply to the Gram- protocol, since Gram+ bacteria contain only a single membrane containing helical integral membrane proteins and not outer membrane beta barrels (OMPbb's) (although it's worth noting that things may be so simple, since some species of Mycobacterium may contain outer membranes and associated OMPs <http://dx.doi.org/10.1016/j.tube.2008.02.004> ). By default, _inmembrane_ uses BOMP (http://services.cbu.uib.no/tools/bomp/) and TMBETADISC-RBF (http://rbf.bioinfo.tw/~sachen/OMPpredict/TMBETADISC-RBF.php).

- `barrel_programs`: a list of web services for outer membrane &beta-barrel prediction ['bomp', 'tmbetadisc-rbf'],
- `bomp_clearly_cutoff`: 3, # if greater than or equal to this, always classify as an OM(barrel)
- `bomp_maybe_cutoff`: 1, # must also have a signal peptide to be OM(barrel)
- `tmbetadisc_rbf_method`: 'aadp', # TMBETADISC-RBF method: aa (single amino acid composition based), dp (dipeptide composition based), aadp (both aa and dp) or pssm (position specific mutation matrix based)

A plugins for TMBETA (`tmbeta`) is available but disabled by default since the general prediction accuracy or web service reliability was found to be lacking during testing. Previously, _inmembrane_ also provided a plugin for TMBHUNT (`tmbhunt`), however this has been removed since the web service is no longer active.

Since topology prediction for arbitrary OMPbb sequences is still very unreliable, the Gram- protocol does not attempt to do loop length prediction but simply annotates potential OMPbb proteins as such.
