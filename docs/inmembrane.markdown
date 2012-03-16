# __inmembrane__, a bioinformatic workflow for annotation of bacterial cell surface proteomes

Andrew J. Perry and Bosco K. Ho
_Department of Biochemistry, Monash University, Melbourne, Australia_

## Abstract 

_inmembrane_ is a tool to predict the surface exposed regions of membrane proteins in sets of bacterial protein sequences. It is intended to be a direct replacement for SurfG+, which originally implemented a protocol for predicting extracellular regions of polypeptide in gram positive bacterial proteomes. In addition, we have implemented a protocol for predicting outer membrane proteomes of gram negative bacteria. _inmembrane_ provides an accessible and transparent code base through use of a modern scripting language. This leads to a program that is much more amenable to modification and extensibility, and provides an example of writing such glue programs for bioinformatic analysis. The program is hosted on a public open-source repository http://github.com/boscoh/inmembrane.

## Introduction

A common task in bioinformatics is to integrate the results of protein prediction programs to deduce complex properties of proteins. In studies of membrane proteomes, quick annotation of an experimentally detected set of the proteins can help highlight sequences of unexpected localization, and can alert researchers to possible contamination from other subcellular fractions. Ultimately, a concise summary of the properties of the detected membrane proteins in a particular proteomic dataset allows meaningful comparisons between different bacterial strains, species, and their responses in membrane remodelling to host and enviromental challenges.

In studies of membrane proteomes, quick annotation of an experimentally detected set of the proteins can help detect sequences of unexpected localization, and can alert researchers to possible contamination from other subcellular fractions. Ultimately, a concise summary of the properties of the detected membrane proteins in a particular proteomic set allows meaningful comparisons between different bacterial strains, species, and their responses in membrane remodelling to host and enviromental challenges.
raw/master/docs/images/membrane_topologies.png
![Membrane topologies](https://github.com/boscoh/inmembrane/raw/master/docs/images/membrane_topologies.png "Figure - Membrane topologies") - topologies of proteins associated with the Gram-negative inner and outer membranes (a), and the Gram-positive plasma membrane and cell wall (b).

>Topologies represented in Gram-negative bacterial inner membrane include (left to right) polytopic transmembrane proteins, monotopic transmembrane proteins and lipoproteins on the periplasmic side of the membrane which are anchored via a lipid moeity covalently attached to the N-terminal cysteine ("CD", where "D" denotes an Asp outer membrane avoidance signal at position 2 (Masuda et al 2002)). The outer membrane also contains lipoproteins, usually on the inner leaflet exposed to the periplasm, however unlike the inner membrane the outer membrane contains beta-barrel membrane proteins ("beta"), some with large extracellular domains exposed on the surface. Akin to the Gram-negative inner membrane, the Gram-positive inner membrane contains mono and polytopic transmembrane proteins and lipoproteins. Gram-positive bacteria also display surface proteins associated covalently or non-covalently with the cell wall peptidoglycan layer via a number of "surface motifs", such as the LPxTG, LysM. Some proteins are also secreted into the extracellular milieu. A subset of Gram-positive bacteria (the Acinetobacterace) have also been shown to contain beta-barrel membrane proteins in their plasma membrane.


One example of this type of annotation is the program SurfG+, which is designed to predict proteins that are exposed on the surface of Gram+ bacteria. SurfG+ is a Java program that carries out batch processing of several standard bioinformatic tools to specifically identify Gram+ bacterial proteins that may be exposed out of the peptidoglycan layer of the bacterium. These predictions are intended to identify a set of proteins that would be amenable to cell-surface protease shaving experiments. SurfG+ itself does not carry out any extensive analysis, but rather relies on a transmembrane helix predictor (_TMMOD_), a secretion signal predictor (_SignalP_), a lipoprotein signal predictior (_LipoP_) and a sequence alignment for protein profiles (_HMMER_). All these programs function as Linux command tools that take FASTA sequences as input and generates output as formatted text.

Nevertheless, _SurfG+_ suffers several problems that plague much bioinformatic software. The most egregrious being that the program is not generally available anymore. The URL mentioned in the original reference no longer exists, even though the paper was only published in 2009. With some searching, we were able to find a source-code repository however due to limited documentation and dependancies which are no longer readily available, we were unable to run _SurfG+_. While Java typically offers the advantage of producing faster executable bytecode than Python (CPython) for equivalent computations, this difference in efficieny is not a major concern for workflows such as _inmembrane_ and _SurfG+_ where the bulk of the CPU intensive operations are carried out by external programs. More critical is the transparency, clarity and modifiability of the code, such that scientists can verify and understand it's behavioir. In contrast to a typical Python program, _SurfG+_, is written in the enterprise Java style of programming, making it far more difficult to read and debug.

Other programs for global prediction of subcellular localization of bacterial proteins exist. Most notable is PSORTb v3.0 (Yu, et al, 2010).

Nevertheless, _SurfG+_ suffers several problems that plague much bioinformatic software. The most egregrious being that the program is not generally available anymore. We were able to find a poorly-managed source-code repository but the we were not able to get the program to work. With the source code, it was possible to deduce the basic algorithm of the program, however, the software architecture of SurfG+ made it extremely difficult to debug properly. Although it may appear at first, to be a redundant effort, there are several compelling reasons for rewriting SurfG+. First, the URL mentioned in the original reference no longer exists, even though the paper was only published in 2009. Second, after a Google search, we managed to find a compiled version of SurfG+ but could not get it to work under Mac OSX or Ubuntu Linux. Third, after more searching, we managed to find a copy of the source-code. Unfortunately, due to the architecture and dependencies of the program, it would have been prohibitively time-consuming to debug the program to work on our systems.

Since the core algorithm in _SurfG+_ is relatively straightforward, we decided to write _inmembrane_ to replicate the functionality of _SurfG+_, but in a modern scripting language. This lead to considerable simplifiction and clarification of the code base. Compared with the _SurfG+_ Java source code of 700K, _inmembrane_ is around 32K of Python code, including additional functionality not offered by _SurfG+_. The smaller and cleaner code case is substantially easier to reuse and repurpose for other users, with a terseness that greatly facilitates modifiability. Here, we discuss the issues involved in writing robust and accessible bioinformatic source code.

## Methodology

### Program Workflow 

The heart of _inmembrane_ is quite simple. It is a wrapper that takes FASTA sequences and sequentially runs a number of sequence analysis programs (HMMER, LipoP, SignalP and TMHMM). These programs share the common feature of taking FASTA sequences as input and generating text output. The bulk of the computation in _inmembrane_ lies in the parsing of the text output. A small amount of analysis is done at the end to integrate the results from the other programs to generate the text output of _inmembrane_. 

As some of the dependent executables of _inmembrane_ are only available on Linux, this unfortunately restricts _inmembrane_ to be fully operational only on the Linux platform. Given the number of dependencies involved, we have provided a comprehensive set of unit tests.

As _inmembrane_ integrates the intermediate output of a large number of programs, there are many potential points of failure. _inmembrane_ saves all intermediate output in a results folder, which allows expert users to diagnose problems at any step in the pipeline and restart the analysis from that point once the problem is resolved.

### Gram-positive protocol

HMMER 3.0 (Robert et al 2011) searches using hidden Markov models (HMM) derived from Pfam and Superfam are used to detect known Gram-positive surface sequence motifs. These include 
LPxTG (Boekhorst et al, 2005) [[PF00746](http://pfam.sanger.ac.uk/family/PF00746) and the HMM used by SurfG+ (Barinov et al 2009)], 
GW repeat domains (Jonquieres et al, 1999) [Superfam models 0040855, 0040856, 0040857], 
peptidoglycan (PG) binding domain (Type 1) (Foster, 1991) [[PF01471](http://pfam.sanger.ac.uk/family/PF01471), [PF08823](http://pfam.sanger.ac.uk/family/PF08823), [PF09374](http://pfam.sanger.ac.uk/family/PF09374)]], 
Choline binding repeats (Janecek et al, 2000), [[PF01473](http://pfam.sanger.ac.uk/family/PF01473)]
LysM domain (Bateman & Bycroft, 2000) [PF01476](http://pfam.sanger.ac.uk/family/PF01476), Cell wall binding domain (Type 2) (Waligora et al, 2001), [[PF04122](http://pfam.sanger.ac.uk/family/PF04122)]
and
S-layer homology domain (Mesnage et al, 2000) [[PF04122](http://pfam.sanger.ac.uk/family/PF04122)]
motifs.

Lipoprotein signals are detected using LipoP (Agnieszka et al 2003), and signal sequences are detected using SignalP (Jannick et al 2004).

The prescence and topology of transmembrane segments in helical membrane proteins is predicted using TMHMM v2.0 (Robel et al, 2005) and/or MEMSAT3 (Jones, 2007). Since MEMSAT3 executes a PSI-BLAST search to gather homologous sequences it is considerably slower than TMHMM, and as such is turned off by default.

_inmembrane_ collates the results of each analysis, and using the predicted topology of the intergral membrane proteins detected predicts potentially surface exposed loops following the algorithm used by SurfG+. By default, external terminal regions longer than 50 residues and external loops longer than 100 residues are considered to be potentially surface exposed. These values were previously experimentally derived based on membrane shaving experiements with S. pyrogenes and may need modification to suit other species with different cell wall thickness (Barinov et al, 2009).

## Discussion

### Public Open Source Repository

Perhaps the single most important step is to distribute our code on an open-source repository Github. We believe that the use of a dedicated repository provides many advantages over the typical strategy of hosting software on an academic server. Github provides excellent code browsing facility, code history, download links, and robust well-defined URL links. These are generally not provided in typical academic releases of software. Github provides excellent usage statistics to measure the impact of the software, which obviates the need for the dreaded login and registration pages.

Most importantly, storing the code in a large, well-supported repository means the source code will remain accessible in the long term, something that historically most academic labs have shown they cannot provide. As well, Github provides excellent facilities to fork a project. That is, if you were to come across an orphaned project that did something useful but had been abandoned a while ago, it is trivial to create a duplicate and make changes.

The algorithms used here are not difficult, and there is no reason to hide the program behind a hybrid academic/commercial license. We have licensed under the BSD license, which liberally allows reuse of the code in any shape or form.

### Scripting Languages 

The virtues of Python as a language for solving problems in life science research have been previously recognized (Bassi, 2007). One potential downside of Python is it's slower execution speed for computationally intensive tasks when compared with compiled languages or Java. Since _inmembrane_ delegates most of the computationally intensive tasks to external programs, the wrapping, text parsing and analysis code in Python does not become a bottleneck in the overall processing speed.

Using a modern scripting language results in much cleaner code, where the advantages arise mostly from the use of standard dynamic language features in Python, which otherwise would require the creation of large complicated objects in Java. Another great advantage of scripting languages is portability, where the source code itself is the executable. The source code can be executed directly without any explicit compilation step. Modification of the source-code directly modifies the program. This dramatically simplifies the process of extending _inmembrane_ for other purposes.

### Simple Data Structures 

While object-oriented programming (OOP) provides advantages when architecting large enterprise systems, it's overuse for small projects can be a disadvantage. In the recommended Enterprise Java style of programming, as used in _SurfG+_, objects are created through several layers of abstract classes where each field in an object needs to be explicitly specified. To change a field in a datastructure, there are at least 6 places in 3 different files where the code that needs to be changed, which severely restricts the ease of modification for those unfamiliar with the code base. Whilst this level of hierarchy is useful in programs that have highly interacting data-structures, this is not the case for _inmembrane_ and would otherwise add unneeded levels of complexity.

Python provides an incredibly powerful standard data structure termed a 'dictionary', which is conceptually similar to a 'hash table' or 'hash map' in other languages. Dictionaries consists of a set of key-value pairs, where keys and values can be any type of data structure - strings, integers, floats, or even other dictionaries. In _inmembrane_, the program data is represented with a flat dictionary called `protein`. Let's say our FASTA file contains the mouse hemeglobin gene with the ID  `'MOUSE_HEME'`. The properties of `MOUSE_HEME` would then be found in `protein['MOUSE_HEME']`, which is itself a dictionary. `protein['MOUSE_HEME']` contains any arbitary number of different properties, also accessed as key-value pairs. For instance, the sequence length of the `'MOUSE_HEME'` sequence would be stored in `protein['MOUSE_HEME']['sequence_length']`. This data structure can capture the results of most potential bioinformatic analyses, where new properties are added to `protein` on the fly. The use of a dynamic flat dictionary avoids much of the boilerplate code involved with an OOP style programming.

### Parsing code is particularly simple

If we use a dictionary to represent our data structure, then the main work in _inmembrane_ of running other programs and processing their text output can be encapsulated into a simple function. 

Using the processing of _SignalP_ as an example, we define a function `run_signalp(protein)` which takes the main protein data structure as input. The function runs SignalP, and then parses the text output. Text processing is very easy to write in Python and the processing for _SignalP_ can be done in about 10 lines of code. As `run_signalp`  cycles through the text output of _SignalP_, then for each protein, say `'MOUSE_HEME'`, if a secretion signal is found, a new property is added: `protein['MOUSE_HEME']['is_signalp'] == True`. 

We can thus abstract the main program loop as running a series of functions of the generic form `run_program(protein)`. This provides a simple API to extend new modules by simply adding new functions that annotates the `protein` dictionary. Depending on the details of how programs are executed, and the format of the output, the functions can get a bit more complicated than `run_signalp`. 

We have exploited other examples of dynamic programming in _inmembrane_ to produce, not only terser code, but clear code. For instance, when using _HMMER_ to match sequence profiles of peptidoglcan binding domains, there is no hard-coding of the sequences in program. Rather, the code dynamically searches the `hmm_profiles` directory for profiles, and iterates the search for each FASTA file across the profiles. This is a particularly robust design, as new profiles can be processed by simply adding them to the directory.

### Web-interface for FASTA file

From time to time, bioinformatics researchers will provide useful sequence analysis tools with a HTML form based front end designed for web browsers, but no official machine readable web API, and no downloadable standalone version of their software. While researchers may neglect to provide these interfaces for a multitude of reasons, for end-users the lack of a standalone version or a web API makes automated use for large scale analyses, such as that carried out by _inmembrane_, somewhat awkward and inconvenient. Several of the published tools for the detection of outer membrane beta-barrel proteins we wished to use as part of the _inmembrane_ workflow only provide a browser based interface, and some only allow submission of a single protein sequence at one time. To solve this problem we chose to implement automated queries to these web intefaces using the _twill_ library (C. Titus Brown, http://twill.idyll.org/), with subsequent parsing of any HTML output using the _BeautifulSoup_ library (Leonard Richardson, http://www.crummy.com/software/BeautifulSoup/).

When writing a wrapper for a new service, commands to interface with the web form can be easily tested directly on the Python commandline, or by using _twill_ itself in interactive mode. This allows for quick prototyping of the scraper, prior to implementation as an _inmembrane_ plugin.

> Figure - commandline example of twill

In it's simplest from, a web service API is essentially an agreement between a service provider and their end-users on a machine readable, predictable and stable interface. Since 'screen scraping' as a method of interfacing with a sequence analysis tool does not use a well defined API with an implicit guarantee of stability, it can be prone to breakage when the format of the HTML form or results page is changed, even slightly. While we believe that the approach taken by _twill_ and the robust parsing provided by _BeautifulSoup_ will prevent many upstream changes breaking these wrappers, inevitably breakage will occur. In this case, the simplicity and ease of modifiability of these wrappers then becomes a key feature that allows expert users to fix them if and when it is required.

### Configuration files and ways of running the program

In many bioinformatic programs written in scripting languages, configuration information is dispersed throughout the program, and users are asked to search through the program and modify the  source code. We believe that this is extremely bad practice as this is highly confusing and encourages a lazy style of programming where the configuration concerns have not been isolated.

Instead, _inmembrane_ reads configuration information from an explicit configuration file `inmembrane.config`, where a default version is auto-generated if it is not initally found. This way the main program can be transferred to different systems without having erroneous configuration information embedded in the source code.

The structure of the configuration file is itself in a Python dictionary form with explicit key-value pairs. Internally, the file can be directly read using the `eval` function directly as a dictionary and passed into the main function `process` that accepts a dictionary of parameters as input. The advantage of taking this approach is that we can provide of running the program that avoids the command-line altogether. We can write a short python program `run_example.py` that defines a new parameters dictionary including the input FASTA file and the output file, then loads _inmembrane_ as a library, and executes the `process` function with these parameters. For an end user, they can simple duplicate the `run_example.py` file, change the parameters, save the file, and assuming that Python is installed on the system, double-click on `run_example.py` to run the program.

## Conclusion

_inmembrane_ provides a clean bioinformatic pipeline to analyze proteomes for proteins that are exposed out of the membrane. It has been written in a  modern style of programming that optimizes readability. It has also been designed to be easily extensible and we sincerely hope that _inmembrane_ will be extended, modified and improved by other researchers. We welcome other researchers to join us on Github.

## References

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

Mesnage, S., Fontaine, T., Mignot, T., Delepierre, M. et al., Bacterial SLH domain proteins are noncovalently anchored to the cell surface via a conserved mechanism involving wall polysaccharide pyruvylation. __EMBO J.__ 2000, 19, 4473–4484. <http://dx.doi.org/10.1093/emboj/19.17.4473>
