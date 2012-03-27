# __inmembrane__, a bioinformatic workflow for annotation of bacterial cell surface proteomes

Andrew J. Perry and Bosco K. Ho
_Department of Biochemistry, Monash University, Melbourne, Australia_

# Abstract 

_inmembrane_ is a tool to predict the surface-exposed regions of membrane proteins in sets of bacterial protein sequences. It is intended to be a direct replacement for SurfG+, which implemented such a protocol for Gram-positive bacterial proteomes. Through the use of a modern scripting language, _inmembrane_ provides a more accessible code base that is easier to modify, and provides a useful example of writing programrams for bioinformatic analysis. The program is hosted on the Github repository http://github.com/boscoh/inmembrane.

# Background

A common task in bioinformatics is to integrate the results of protein prediction programs to deduce complex properties of proteins. In studies of membrane proteomes, quick annotation of an experimentally detected set of the proteins can help highlight sequences of unexpected localization, and can alert researchers to possible contamination from other subcellular fractions. Ultimately, a concise summary of the properties of the detected membrane proteins in a particular proteomic dataset allows meaningful comparisons between different bacterial strains, species, and their responses in membrane remodelling to host and enviromental challenges.

In studies of membrane proteomes, quick annotation of an experimentally detected set of the proteins can help detect sequences of unexpected localization, and can alert researchers to possible contamination from other subcellular fractions. Ultimately, a concise summary of the properties of the detected membrane proteins in a particular proteomic set allows meaningful comparisons between different bacterial strains, species, and their responses in membrane remodelling to host and enviromental challenges.

*****
![Membrane topologies](https://github.com/boscoh/inmembrane/raw/master/docs/images/membrane_topologies.png "Figure - Membrane topologies")

>Topologies represented in Gram-negative bacterial inner membrane include (left to right) polytopic transmembrane proteins, monotopic transmembrane proteins and lipoproteins on the periplasmic side of the membrane which are anchored via a lipid moeity covalently attached to the N-terminal cysteine ("CD", where "D" denotes an Asp outer membrane avoidance signal at position 2 (Masuda et al 2002)). The outer membrane also contains lipoproteins, usually on the inner leaflet exposed to the periplasm, however unlike the inner membrane the outer membrane contains beta-barrel membrane proteins ("beta"), some with large extracellular domains exposed on the surface. Akin to the Gram-negative inner membrane, the Gram-positive inner membrane contains mono and polytopic transmembrane proteins and lipoproteins. Gram-positive bacteria also display surface proteins associated covalently or non-covalently with the cell wall peptidoglycan layer via a number of "surface motifs", such as the LPxTG, LysM. Some proteins are also secreted into the extracellular milieu. A subset of Gram-positive bacteria (the Acinetobacterace) have also been shown to contain beta-barrel membrane proteins in their plasma membrane.
*****

A number of published software packages exist for global prediction of subcellular localization of bacterial proteins exist. Most notable is _PSORTb v3.0_ (Yu, et al, 2010) which predicts general subcellular localization for Gram-positive, Gram-negative and Archaeal proteins sequences. _Augur_ (__ref__)is a specialized tool which enables prediction and proteome comparison of Gram-positive surface proteins via a web interface. LocateP (__ref__) is a pipeline wrapping existing specific localization predictors, and provides a web accessible database of pre-calculated subcellular localization for Gram-positive proteomes. The publicly available LocateP server does not allow analysis of arbitrary sets of sequences. While the source code for both _PSORTb 3.0_ is available under an open source license, the code for the other packages is not generally available for download.

An extension to general membrane localization prediction is the analysis of membrane protein topology to identify those with prominent surface exposed loops. These potentially surface exposed (PSE) proteins constitute attractive vaccine candidates. One such workflow for annotation of PSE proteins is the program SurfG+, which focuses on Gram-positive bacterial proteomes. SurfG+ is a Java program that carries out batch processing of several standard bioinformatic tools to specifically predict proteins that  protude out of the peptidoglycan layer of the bacterium. These predictions are intended to identify a set of proteins that would be amenable to cell-surface protease shaving experiments. SurfG+ itself does not carry out any extensive analysis itself, but rather leverages the results of a transmembrane helix predictor (_TMMOD_) (Robel et al, 2005), a secretion signal predictor (_SignalP_) (Jannick et al 2004), a lipoprotein signal predictior (_LipoP_) (Agnieszka et al 2003) and a sequence alignment for protein profiles (_HMMER_) (Robert et al 2011). 

Nevertheless, _SurfG+_ suffers several problems that plague much bioinformatic software. Despite being published in 2009, the URL mentioned in the original reference no longer exists. We were able to find a [source-code repository](https://mulcyber.toulouse.inra.fr/projects/surfgplus) but the we were not able to get the program to work, due in part to dependencies that are not longer generally available for download.

Since the core algorithm in _SurfG+_ is relatively straightforward, we decided to replicate the functionality of _SurfG+_ by writing _inmembrane_ in a modern scripting language. This lead to considerable simplifiction and clarification of the code base. Compared with the _SurfG+_ Java source code of 700K, _inmembrane_, without dependencies, is around 70K of Python code and includes additional functionality not offered by _SurfG+_. The smaller code case is substantially easier to reuse and repurpose for other users. Here, we discuss the issues involved in writing robust and accessible bioinformatic source code.

## Methods and Implementation
_inmembrane_ is primarily designed to be run locally via the command line. The input is a set of sequences in FASTA format, the output is plain text, including a summary table and comma-separated-value (CSV) format suitable for import into spreadsheet software or scripted text processing.

> __Figure__ - snippet of example output ?

A set of unit tests, executable via the _run_test.py_ script enables users and developers to quickly verify if their _inmembrane_ installation, with dependencies, is functioning as expected.

### Gram-positive protocol
The _inmembrane_ Gram-positive surface protocol leverages a number of existing single localization predictors, including transmembrane topology prediction, to deduce a likely subcellular localization and expected surface exposure of a given proteome. Each sequence is annotated by every predictor, and then these annotations are used by the business logic of _inmembrane_ to classify proteins as potentially surface exposed ("PSE"), "Secreted", or the non-exposed classes "Cytoplasmic" and "Membrane".

Annotations applied are as follows. HMMER 3.0 (Robert et al 2011) searches using hidden Markov models (HMM) derived from Pfam and Superfam are used to detect known Gram-positive surface sequence motifs. These include 
LPxTG (Boekhorst et al, 2005) [[PF00746](http://pfam.sanger.ac.uk/family/PF00746) and the HMM used by SurfG+ (Barinov et al 2009)], 
GW repeat domains (Jonquieres et al, 1999) [Superfam models 0040855, 0040856, 0040857], 
peptidoglycan (PG) binding domain (Type 1) (Foster, 1991) [[PF01471](http://pfam.sanger.ac.uk/family/PF01471), [PF08823](http://pfam.sanger.ac.uk/family/PF08823), [PF09374](http://pfam.sanger.ac.uk/family/PF09374)]], 
Choline binding repeats (Janecek et al, 2000), [[PF01473](http://pfam.sanger.ac.uk/family/PF01473)]
LysM domain (Bateman & Bycroft, 2000) [PF01476](http://pfam.sanger.ac.uk/family/PF01476), Cell wall binding domain (Type 2) (Waligora et al, 2001), [[PF04122](http://pfam.sanger.ac.uk/family/PF04122)],
S-layer homology domain (Mesnage et al, 2000) [[PF04122](http://pfam.sanger.ac.uk/family/PF04122)]
motifs and the NLPC_P60 cell wall associated domain (Anantharaman & Aravind, 2003) [[PF00877](http://pfam.sanger.ac.uk/family/PF00877)]. PFAM HMMs are from most recent version of at the time of writing, release 26.0.

Lipoprotein signals are detected using LipoP (Agnieszka et al 2003), and signal sequences are detected using SignalP (Jannick et al 2004), including detection of signal peptidase cleavage sites.

The prescence and topology of transmembrane segments in helical membrane proteins is predicted using TMHMM v2.0 (Robel et al, 2005) and/or MEMSAT3 (Jones, 2007). Since MEMSAT3 executes a PSI-BLAST search to gather homologous sequences it is considerably slower than TMHMM, and as such is turned off by default.

_inmembrane_ collates the results of each analysis, and using the predicted topology of the intergral membrane proteins detected, predicts potentially surface-exposed loops following the algorithm used by SurfG+. By default, external terminal regions longer than 50 residues and external loops longer than 100 residues are considered to be potentially surface exposed. These values were previously experimentally derived based on membrane shaving experiements with _S. pyrogenes_ and may need modification to suit other species with different cell wall thickness (Barinov et al, 2009).

### Gram-negative protocol

>__TODO__

### Future protocols

__inmembrane__ is designed such that new workflows for annotation of membrane proteomes can be added relatively easily. Wrappers for programs that annotate a sequence with a particular feature can be added to `inmembrane/plugins/` following the example of existing plugins. The `inmembrane/plugin/signalp4.py` and `inmembrane/plugin/lipop1.py` plugins provide good templates for adoption and modification. In the simplest case, this means that if a superior method for signal peptide, transmembrane segment or lipoprotein prediction is developed, it will be straightforward to write a new plugin wrapping it for inclusion in the protocol, either as a parallel analysis or replacing one of the existing predictors. New _protocols_ can be added to the `inmembrane/protocols` directory, and selected for execution by changing _protocol_ parameter in the `inmembrane.config` file. Currently, we have implemented two protocols, _gram\_pos_, for prediction of PSE proteins in Gram-positive bacteria, and _gram\_neg_, for general annotation of Gram-positive subcellular localization.

## Discussion

### Public Open Source Repository
Perhaps the single most important step in improving code is to distribute it on an open-source repository such as Github. We believe that the use of a dedicated repository provides many advantages over the typical strategy of hosting software on an academic server. Github provides excellent code-browsing facility, code history, download links, and robust well-defined URL links. These are generally not provided in typical academic releases of software. Github provides excellent usage statistics to measure the impact of the software, which obviates the need for the dreaded login and registration pages.

Most importantly, storing the code in a well-supported repository means the source code will remain accessible in the long term, something that historically most academic labs have shown they cannot provide. If you were to come across an abandonned project on Github, it is trivial to fork the project - create a duplicate project to change and/or improve the program. Since the algorithms used here are rather straightforward, we have used a liberal BSD license, rather than the usual hybrid academic/commercial license, which allows reuse of the code in any form.

### Program setup and workflow

The heart of _inmembrane_ is simple: it takes FASTA sequences, sequentially runs the sequences through a number of sequence analysis programs (HMMER, LipoP, SignalP, TMHMM and optionally MEMSAT3), and generates the results in a text output. The bulk of the computation in _inmembrane_ lies in the parsing of the text output of the external programs where a small amount of analysis is done at the end.

As _inmembrane_ integrates the  output of a large number of external dependencies, there are many potential points of failure. As such, _inmembrane_ saves all intermediate output into a results folder, and a comprehensive set of unit tests is provided to help diagnose dependencies issues. Unfortunately, as some of the dependent executables of _inmembrane_ are only available on Linux, the full operation of _inmembrane_ is restricted to the Linux platform.

In many bioinformatic programs written in scripting languages, configuration information is dispersed throughout the program, and users are asked to search through the program and modify the  source code. We believe this is extremely bad practice as this is highly confusing and encourages a lazy style of programming where the configuration concerns have not been isolated. Instead, _inmembrane_ reads configuration information from an explicit configuration file `inmembrane.config`, where a default version is auto-generated if it is not initally found.

We have structured _inmembrane_ such that you can write a short Python script that incorporate a specific configuration and executes _inmembrane_ directly. This provides a convenient record of each individual analysis, as well as file that can be executed through a file-manager by double-clicking.

### Scripting Languages 

The virtues of Python as a language for solving problems in life science research have been previously recognized (Bassi, 2007). One potential downside of Python is it's slower execution speed for computationally intensive tasks when compared with compiled languages or Java. Since _inmembrane_ delegates most of the computationally intensive tasks to external programs, the wrapping, text parsing and analysis code in Python does not become a bottleneck in the overall processing speed.

Programs written in Java often follow an object-oriented programming (OOP) approach. Although OOP provides advantages when architecting large enterprise systems, it's overuse for small projects can be a disadvantage. In the recommended Enterprise Java style of programming used in _SurfG+_, objects are created through several layers of abstract classes where each field in an object needs to be explicitly specified. To change a field in a data structure, there are at least 6 places in 3 different files where the code that needs to be changed, which severely restricts the ease of modification for those unfamiliar with the code base. Whilst this level of hierarchy is useful in programs that have highly interacting data-structures, this is not the case for the bioinformatic analysis here and adds otherwise unneeded levels of complexity.

Using a modern scripting language such as Python results in  cleaner code, where the advantages arise mostly from the use of standard dynamic language features, which otherwise would require the creation of large complicated objects in Java. Another advantage is portability, where the Python source code itself is directly executable. This avoids the often complicated steps involved in compilation when modifiying the source-code. 

### Simple Data Structures Allows Simple Text-Parsing

In _inmembrane_, the standard Python dictionary is used that provides a flexible way to represent data and allows extremely simple parsing code to be written. The Python 'dictionary', which is conceptually similar to a 'hash table' or 'hash map' in other languages, consists of a set of key-value pairs, where keys and values can be any type of data structure - strings, integers, floats, or even other dictionaries. 

The main program data in _inmembrane_ is represented with a flat dictionary called `protein`. Let's say our FASTA file contains the mouse hemeglobin gene with the ID  `'MOUSE_HEME'`. The properties of `MOUSE_HEME` would then be found in `protein['MOUSE_HEME']`, which is itself a dictionary. `protein['MOUSE_HEME']` contains any arbitary number of different properties, also accessed as key-value pairs. For instance, the sequence length of the `'MOUSE_HEME'` sequence would be stored in `protein['MOUSE_HEME']['sequence_length']`. This data structure can capture the results of most potential bioinformatic analyses, where new properties are added to `protein` on the fly. The use of a dynamic flat dictionary avoids much of the boilerplate code involved with an OOP style programming.

If we use a dictionary to represent our data structure, then the main work in _inmembrane_ of running other programs and processing their text output can be encapsulated into a simple function. For example with _SignalP_, we define a function `annotate_signalp(params, protein)` which takes the main protein data structure as input. The function runs SignalP, and then parses the text output. Text processing is very easy to write in Python and the processing for _SignalP_ can be done in about 10 lines of code. As `annotate_signalp` cycles through the text output of _SignalP_, then for each protein, if a secretion signal is found, a new property is added: `protein['MOUSE_HEME']['is_signalp'] == True`. 

We can thus abstract the main program loop as running a series of functions of the generic form `annotate_program(params, protein)`. This provides a simple API to extend new modules by simply adding new functions that annotates the `protein` dictionary. 

Another example of cleaner code through dynamic programming is in the _HMMER_ peptide motif matching. Instead of hard-coding the sequence profiles search as in _SurfG+_, _inmembrane_ dynamically searches the `hmm_profiles_dir` directory for profiles to search against. Conveniently, new profiles can then be processed by simply dropping them into this directory.

### Tests with Gram-positive bacteria

The field of bioinformatics changes quickly, and in the few years between the release of SurfG+, some of the software used in SurfG+ is no longer available. As a result we could not use exactly the same versions of the binaries used in SurfG+. For instance TMMOD is no longer released as a binary and SignalP has progressed to Version 4.0. Nevertheless, _inmembrane_ produces comparable results to SurfG+ for the 4 bacterial genomes:

<pre>
             S.pyogenes L.acidophilus L.johnsonii   L.gasseri  L.bulgaricus
Accession    AE004092    CP000033      AE017198      CP000413    CR954253
Program       S    i       S   i         S   i        S   i       S   i
Cytoplasmic  1243 1234   1290 1280     1248 1234    1262 1240   1132 1119                         
Membrane      236 239     315 332       357 358      298 303     244 264                          
PSE           140 176     169 189       176 204      157 189     116 138                           
Secreted       78 47       88 61         40 25        38 23       70 41                    
Total        1696 1696   1862 1862     1821 1821    1755 1755   1562 1562
</pre>

Columns labelled 'S' are _SurfG+_ results and 'i' are _inmembrane_ results.

This can be compared to PSortB 3.0 classifcation of the organisms.

<pre>
               S.pyogenes L.acidophilus L.johnsonii L.gasseri L.bulgaricus
Cellwall            24          46           26          27         19
Cytoplasmic         884         855          826         804        743
Cytopl. Membrane    432         519          548         489        440
Extracellular       28          32           16          13         15
Unknown             323         402          394         419        307
Unknown/multiple    5           8            11          3          5
Total               1696        1862         1821        1755       1529
</pre>

### Interfacing with "web services" for further analysis

From time to time, bioinformatics researchers will provide useful sequence analysis tools with a HTML form based front end designed for web browsers, but no official machine readable web API, and no downloadable standalone version of their software. While researchers may neglect to provide these interfaces for a multitude of reasons, for end-users the lack of a standalone version or a web API makes automated use for large scale analyses, such as that carried out by _inmembrane_, somewhat awkward and inconvenient. Several of the published tools for the detection of outer membrane beta-barrel proteins we wished to use as part of the _inmembrane_ workflow only provide a browser based interface, and some only allow submission of a single protein sequence at one time. To solve this problem we chose to implement automated queries to these web intefaces using the _twill_ library (C. Titus Brown, http://twill.idyll.org/), with subsequent parsing of any HTML output using the _BeautifulSoup_ library (Leonard Richardson, http://www.crummy.com/software/BeautifulSoup/).

When writing a wrapper for a new service, commands to interface with the web form can be easily tested directly on the Python commandline, or by using _twill_ itself in interactive mode. This allows for quick prototyping of the scraper, prior to implementation as an _inmembrane_ plugin.

> Figure - commandline example of twill

In it's simplest from, a web service API is essentially an agreement between a service provider and their end-users on a machine readable, predictable and stable interface. Since 'screen scraping' as a method of interfacing with a sequence analysis tool does not use a well defined API with an implicit guarantee of stability, it can be prone to breakage when the format of the HTML form or results page is changed, even slightly. While we believe that the approach taken by _twill_ and the robust parsing provided by _BeautifulSoup_ will prevent many upstream changes breaking these wrappers, inevitably breakage will occur. In this case, the simplicity and ease of modifiability of these wrappers then becomes a key feature that allows expert users to fix them if and when it is required.

# Conclusion

_inmembrane_ provides a clean bioinformatic pipeline to analyze proteomes for proteins that are exposed out of the membrane. It has been written in a  modern style of programming that optimizes readability. It has also been designed to be easily extensible and we sincerely hope that _inmembrane_ will be modified and improved by other researchers. We welcome other researchers to join us on Github.

# Availability and requirements

__Project name:__ inmembrane

__Project home page:__ http://github.com/boscoho/inmembrane

__Operating systems:__ Linux

__Programming language:__ Python

__Other requirements:__ HMMER, SignalP, LipoP, TMHMM or MEMSAT3. An Internet connection is required for web services such as BOMP and TMBETA-NET.

__Licence:__ BSD Licence (2-clause)

Any restrictions to use by non-academics: Use of _inmembrane_ itself is unrestricted, however many of the dependencies require special licencing for non-academic use.


# References

Barinov A, Loux V, Hammani A, Nicolas P, Langella P, et al. (2009) Prediction of surface exposed proteins in Streptococcus pyogenes, with a potential application to other Gram-positive bacteria. __Proteomics__ 9: 61-73. <http://dx.doi.org/10.1002/pmic.200800195>

Bassi S (2007) A Primer on Python for Life Science Researchers. __PLoS Comput Biol__ 3(11): e199. <http://dx.doi.org/10.1371/journal.pcbi.0030199>

Jannick Dyrløv Bendtsen, Henrik Nielsen, Gunnar von Heijne and Søren Brunak. (2004) Improved Prediction of Signal Peptides: SignalP 3.0. __J. Mol. Biol.__ 340:783–795.

Robert D. Finn, Jody Clements and Sean R. Eddy. (2011) HMMER web server: interactive sequence similarity searching. __Nucleic Acids Research__ 39:W29–W37.

Agnieszka S. Juncker, Hanni Willenbrock, Gunnar Von Heijne, Søren Brunak, Henrik Nielsen, And Anders Krogh. (2003) Prediction of lipoprotein signal peptides in Gram-negative bacteria. __Protein Science__ 12:1652–1662.

Anders Krogh, Björn Larsson, Gunnar von Heijne and Erik L. L. Sonnhammer (2001) Predicting Transmembrane Protein Topology with a Hidden Markov Model: Application to Complete Genomes. __J. Mol. Biol.__ 305:567-580.

Robel Y. Kahsay1, Guang Gao1 and Li Liao1. An improved hidden Markov model for transmembrane protein detection and topology prediction and its applications to complete genomes (2005) __Bioinformatics__ 21: 1853-1858.

﻿Jones DT. Improving the accuracy of transmembrane protein topology prediction using evolutionary information. __Bioinformatics (Oxford, England)__ 2007 Mar;23(5):538-44. <http://dx.doi.org/10.1093/bioinformatics/btl677>

N.Y. Yu, J.R. Wagner, M.R. Laird, G. Melli, S. Rey, R. Lo, P. Dao, S.C. Sahinalp, M. Ester, L.J. Foster, F.S.L. Brinkman (2010) PSORTb 3.0: Improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes, __Bioinformatics__ 26(13):1608-1615 <http://dx.doi.org/10.1093/bioinformatics/btq249>

﻿Masuda K, Matsuyama S-ichi, Tokuda H. Elucidation of the function of lipoprotein-sorting signals that determine membrane localization. __Proceedings of the National Academy of Sciences of the United States of America__ 2002 May;99(11):7390-5. <http://dx.doi.org/10.1073/pnas.112085599>

Boekhorst, J., de Been, M. W., Kleerebezem, M., Siezen, R. J., Genome-wide detection and analysis of cell wall-bound proteins with LPxTG-like sorting motifs. __J. Bacteriol.__ 2005, 187, 4928–4934. <http://dx.doi.org/10.1128/​JB.187.14.4928-4934.2005>

Jonquieres, R., Bierne, H., Fiedler, F., Gounon, P., Cossart, P., Interaction between the protein InlB of Listeria mono- cytogenes and lipoteichoic acid: A novel mechanism of
protein association at the surface of Gram-positive bacteria. __Mol.Microbiol.__ 1999, 34, 902–914. <http://dx.doi.org/10.1046/j.1365-2958.1999.01652.x>

Foster, S. J., Cloning, expression, sequence analysis and biochemical characterization of an autolytic amidase of Bacillus subtilis 168 trpC2. __J. Gen. Microbiol.__ 1991, 137, 1987–1998. <http://dx.doi.org/10.1099/00221287-137-8-1987>

Janecek, S., Svensson, B., Russell, R. R., Location of repeat elements in glucansucrases of Leuconostoc and Strepto- coccus species. __FEMS Microbiol. Lett.__ 2000, 192, 53–57. <http://dx.doi.org/10.1111/j.1574-6968.2000.tb09358.x>

Bateman, A., Bycroft, M., The structure of a LysM domain from E. coli membrane-bound lytic murein transglycosylase D (MltD). __J. Mol. Biol.__ 2000, 299, 1113–1119. <http://dx.doi.org/10.1006/jmbi.2000.3778>

Waligora, A. J., Hennequin, C., Mullany, P., Bourlioux, P. et al., Characterization of a cell surface protein of Clostridium difficile with adhesive properties. __Infect. Immun.__ 2001, 69, 2144–2153. <http://dx.doi.org/10.1128/​IAI.69.4.2144-2153.2001>

N.Y. Yu, J.R. Wagner, M.R. Laird, G. Melli, S. Rey, R. Lo, P. Dao, S.C. Sahinalp, M. Ester, L.J. Foster, F.S.L. Brinkman (2010) PSORTb 3.0: Improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes, __Bioinformatics__ 26(13):1608-1615 <http://dx.doi.org/10.1093/bioinformatics/btq249>

Mesnage, S., Fontaine, T., Mignot, T., Delepierre, M. et al., Bacterial SLH domain proteins are noncovalently anchored to the cell surface via a conserved mechanism involving wall polysaccharide pyruvylation. __EMBO J.__ 2000, 19, 4473–4484. <http://dx.doi.org/10.1093/emboj/19.17.4473>
