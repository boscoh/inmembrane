# Customization of _inmembrane_ parameters

All adjustable parameters of the running of _inmembrane_ have been collected in a parameters file, which is typically stored in the `inmembrane.config`. On the first use of the program, when no such `inmembrane.config` is found, a default `inmembrane.config` will be generated, otherwise existing `inmembrane.config` will be used.

The default parameters are as following, they will be explained in sections below:

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
    # 'barrel_programs': ['bomp', 'tmbetadisc-rbf', 'tmbhunt', 'tmbeta'],
      'bomp_clearly_cutoff': 3, # >= to this, always classify as an OM(barrel)
      'bomp_maybe_cutoff': 1, # must also have a signal peptide to be OM(barrel)
      'tmbhunt_clearly_cutoff': 0.95,
      'tmbhunt_maybe_cutoff': 0.5,
      'tmbetadisc_rbf_method': 'aadp', # aa, dp, aadp or pssm
    }

Existing `inmembrane.config` files can be changed. The parameters will be explained below, grouped in terms of related parameters.

### General Parameters

These are the general parameters to the program, representing the top-down organization of the work-flow. 

A fasta program is the input, following one of currently two different protocols - one for Gram+ bacteria, and one for Gram- bacteria. The results are saved to a csv file as an output.

- 'fasta': the pathname of the fasta file holding the sequences
- 'protocol': switches between different protocols for analysising the proteomes
- 'csv': the name of the output file in CSV (Excel-compatible) format. Defaults to a file using the same filename as 'fasta'

For potential debugging, it's important that all intermediate files are saved, and they are saved in: 

- 'out_dir': the directory where all intermediate output files are located

### Optional - Binaries or Web-based API's

One of the problems that _inmembrane_ tackles is the complexity of using and installing multiple binaries. In older programs, the only option was to install local binaires of the external programs. This raises problems involved with getting the binary for the right platform, compilation with the right tools, and installation in the correct file system. In contrast, many of these tools now provide a web-interface. In _inmembrane_, we have added facilities to use the web service, whether through a RESTful interface, or through a SOAP interface, when provided.

Consequently, the parameters for the binaries are needed, where _inmembrane_ will attempt to run the program in the default directory, or if a full-path name is given, _inmembrane_ can run the binary in the correct location. 
In order to turn off the local binary option, if in a *_bin parameter is given the empty string, the local binary step will be skipped and the web version will be used instead.

- 'signalp4_bin', 'lipop1_bin', 'tmhmm_bin', 'memsat3_bin': the full pathname of the binaries used for the analysis. If empty, it is not run, and the web version is used instead

### The Surface-Exposure Algorithm

The key part of the surface-exposure prediction evaluates the length of loops in the extracellular part of membrane proteins. This evaluation requires the key parameter of loop cutoff that determines the minimum length of a loop that is required to potentially pass through the glycan layer. 

In the original _SurfG+_ paper, for Gram+ bacteria, it was found that a loop lenght of ~50 amino acids is required. Therefore for internal loops, this gives a length ~100 amino acids, which allows the loop to pass through the layer and back down again.

- 'terminal_exposed_loop_min': this is the length in terms of amino acids given to define an external loop of a protein to be long enough to stick out of the bacterial cell-wall to be considered surface-epxosed
- 'internal_exposed_loop_min': similarly to above, except is double the length as internal loops need to go up, and down again, for continuity with the rest of the protein.

There are several different way to identify proteins that have loops that potentially stick out of the glycan cell-wall of the bacteria. These will be discussed in the following sections.

### Proteins attached to the Glycan Cell-wall and the Lipid bilayer

The clearest example of proteins that may stick of the glycan layer are proteins that possess known cell-wall attachment motifs. 

Motifs for motifs that bind to the glycan cell-wall are taken from PFAM as HMMER profiles, and stored into a directory included in the _inmembrane_ package. In the default `inmembrane.config` file, the location of this directory will be put in `hmm_profiles_dir`. 

Subsequently, `HMMER` will be used to search the profiles against the input query sequences. The cutoffs for the identification of a motif is given by `hmm_evalue_max` and `hmm_score_min`.

- 'hmmsearch3_bin': the binary for HMMER
- 'hmm_profiles_dir': the directory that holds the HMMER profiles for cell-wall signals
- 'hmm_evalue_max': the maximum expectation value cutoff for acceptable prediction in HMMER
- 'hmm_score_min': the minimum score cutoff for acceptatble prediction in HMMER

Another source of proteins attached to the membrane are lipo-proteins. LipoP is used to identify motifs for the cell to attach a protein to the lipid layer of the protein. 

- 'lipop_bin': the location of the binary of LipoP, or if this is given as '', then it uses the web-version at http://http://www.cbs.dtu.dk/services/LipoP/.

### Transmembrane-&alpha;-helix Predictors

The other main source of proteins in Gram+ proteins with external loops are the identification of proteins with transmembrane &alpha-helices. The internal loops of proteins with multiple transmembrane helices, as well as the terminal loops are all candidates to be exposed out of the glycan layer.

There is no dominant transmembrane helix predictor, and thus _inmembrane_ has been written to allow pluggable options for transmembrane helix prediction. To allow for this, _inmembrane_ allows the results of more than one helix predictor, and the one that give the longest loops that extends out of the extracellular region will be used. As such the parameter 'helix_programs' lists all the transmembrane helix predictors to be used. Currently, 

- 'helix_programs': a list of transmembrane-helix predictors to be used. Currently, _inmembrane_ can handle `tmhmm` and `memsat`.

### &beta;-barrel predictors

- 'barrel_programs': list of programs to do &beta-barel ['bomp', 'tmbetadisc-rbf'],
- 'bomp_clearly_cutoff': 3, # >= to this, always classify as an OM(barrel)
- 'bomp_maybe_cutoff': 1, # must also have a signal peptide to be OM(barrel)
- 'tmbhunt_clearly_cutoff': 0.95,
- 'tmbhunt_maybe_cutoff': 0.5,
- 'tmbetadisc_rbf_method': 'aadp', # aa, dp, aadp or pssm


