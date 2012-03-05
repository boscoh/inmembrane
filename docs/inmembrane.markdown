# __inmembrane__, a bioinformatic workflow for annotation of bacterial cell surface proteomes

Andrew J. Perry and Bosco K. Ho  
_Department of Biochemistry, Monash University, Melbourne, Australia_

## Abstract 

_inmembrane_ is a tool to predict the surface-exposed regions of membrane proteins in sets of bacterial protein sequences. It is intended to be a direct replacement for SurfG+, which implemented such a protocol for Gram+ bacterial proteomes. Through the use of a modern scripting language, _inmembrane_ provides a more accessible code base that is easier to modify, and provides a useful example of writing programs for bioinformatic analysis. The program is hosted on the github repository http://github.com/boscoh/inmembrane.

## Introduction

A common task in bioinformatics is to integrate the results of protein prediction programs to deduce complex properties of proteins. As sequencing techniques improve, the burgeoning volume of data will require a greater capacity to analyze such sequences in an automated way.

_inmembrane_ is a tool to predict the surface exposed regions of membrane proteins in sets of bacterial protein sequences. It is intended to be a direct replacement for SurfG+, which originally implemented a protocol for predicting extracellular regions of polypeptide in gram positive bacterial proteomes. In addition, we have implemented a protocol for predicting outer membrane proteomes of gram negative bacteria. _inmembrane_ provides an accessible and transparent code base through use of a modern scripting language. This leads to a program that is much more amenable to modification and extensibility, and provides an example of writing such glue programs for bioinformatic analysis. The program is hosted on a public open-source repository http://github.com/boscoh/inmembrane.

## Introduction

A common task in bioinformatics is to integrate the results of protein prediction programs to deduce complex properties of proteins. In studies of membrane proteomes, quick annotation of an experimentally detected set of the proteins can help highlight sequences of unexpected localization, and can alert researchers to possible contamination from other subcellular fractions. Ultimately, a concise summary of the properties of the detected membrane proteins in a particular proteomic dataset allows meaningful comparisons between different bacterial strains, species, and their responses in membrane remodelling to host and enviromental challenges.

In studies of membrane proteomes, quick annotation of an experimentally detected set of the proteins can help detect sequences of unexpected localization, and can alert researchers to possible contamination from other subcellular fractions. Ultimately, a concise summary of the properties of the detected membrane proteins in a particular proteomic set allows meaningful comparisons between different bacterial strains, species, and their responses in membrane remodelling to host and enviromental challenges.

One example of this type of annotation is the program SurfG+, which is designed to predict proteins that are exposed on the surface of Gram+ bacteria. SurfG+ is a Java program that carries out batch processing of several standard bioinformatic tools to specifically identify Gram+ bacterial proteins that may be exposed out of the peptidoglycan layer of the bacterium. These predictions are intended to identify a set of proteins that would be amenable to cell-surface protease shaving experiments. SurfG+ itself does not carry out any extensive analysis, but rather relies on a transmembrane helix predictor (_TMMOD_), a secretion signal predictor (_SignalP_), a lipoprotein signal predictior (_LipoP_) and a sequence alignment for protein profiles (_HMMER_). All these programs function as Linux command tools that take FASTA sequences as input and generates output as formatted text.

Nevertheless, _SurfG+_ suffers several problems that plague much bioinformatic software. The most egregrious being that the program is not generally available anymore. The URL mentioned in the original reference no longer exists, even though the paper was only published in 2009. We were able to find a poorly-managed source-code repository but we could not get the program to work. The program was written in the enterprise Java style of programming, making it difficult to read, and almost impossible to debug. 

Other programs for global prediction of subcellular localization of bacterial proteins exist. Most notable is PSORTb v3.0 (Yu, et al, 2010).

Nevertheless, _SurfG+_ suffers several problems that plague much bioinformatic software. The most egregrious being that the program is not generally available anymore. We were able to find a poorly-managed source-code repository but the we were not able to get the program to work. With the source code, it was possible to deduce the basic algorithm of the program, however, the software architecture of SurfG+ made it extremely difficult to debug properly. Although it may appear at first, to be a redundant effort, there are several compelling reasons for rewriting SurfG+. First, the URL mentioned in the original reference no longer exists, even though the paper was only published in 2009. Second, after a Google search, we managed to find a compiled version of SurfG+ but could not get it to work under Mac OSX or Ubuntu Linux. Third, after more searching, we managed to find a copy of the source-code. Unfortunately, due to the architecture and dependencies of the program, it would have been prohibitively time-consuming to debug the program to work on our systems.

However, at it's core _SurfG+_ uses a relatively straightforward algorithm and consequently, we decided to write a tool to reproduce the functionality of _SurfG+_. In the process of improving the architecture of SurfG+ and writing _inmembrane_, we have generated a better template for this class of bioinformatic software. One example is that from a Java source of 700K, we have reduced it down ~32K of Python code, where the terseness of the code greatly facilitates modifiability. Here, we discuss the issues involved in writing robust and comprehensible bioinformatic source code.

Since the core algorithm in _SurfG+_ is relatively straightforward, we decided to write _inmembrane_ to replicate the functionality of _SurfG+_, but in a modern scripting language. This lead to a major simplifiction of the code base. As a rough measure of the simplification, we reduced the original Java source of 700K to 10K of Python code. Although none of the techniques used in _inmembrane_ is particularly new, they may be unfamiliar to the bioinformatic community. Here, we discuss the ideas that lead to a smaller and cleaner code case that is substantially easier to reuse and repurpose for other users.

## Discussion

### Public Open Source Repository

Perhaps the single most important step is to distribute our code on an open-source repository Github. We believe that the use of a dedicated repository provides many advantages over the typical strategy of hosting software on an academic server. Github provides excellent code browsing facility, code history, download links, and robust well-defined URL links. These are generally not provided in typical academic releases of software. Github provides excellent usage statistics to measure the impact of the software, which obviates the need for the dreaded login and registration pages.

Most importantly, storing the code in a large, well-supported repository means the source code will remain accessible in the long term, something that historically most academic labs have shown they cannot provide. As well, Github provides excellent facilities to fork a project. That is, if you were to come across an orphaned project that did something useful but had been abandoned a while ago, it is trivial to create a duplicate and make changes.

The algorithms used here are not difficult, and there is no reason to hide the program behind a hybrid academic/commercial license. We have licensed under the BSD license, which liberally allows reuse of the code in any shape or form.

### Program Workflow 

The heart of _inmembrane_ is quite simple. It is a wrapper that takes FASTA sequences and sequentially runs a number of sequence analysis programs (HMMER, LipoP, SignalP and TMMOD). These programs share the common feature of taking FASTA sequences as input and generating text output. The bulk of the computation in _inmembrane_ lies in the parsing of the text output. A small amount of analysis is done at the end to integrate the results from the other programs to generate the text output of _inmembrane_. 

As some of the dependent executables of _inmembrane_ are only available on Linux, this unfortunately restricts _inmembrane_ to be fully operational only on the Linux platform. Given the number of dependencies involved, we have provided a comprehensive set of unit tests.

As _inmembrane_ integrates the intermediate output of a large number of programs, there are many potential points of failure. _inmembrane_ saves all intermediate output in a results folder, which allows expert users to diagnose problems with the dependencies.

### Scripting Languages 

The virtues of Python as a language for solving problems in life science research have been previously recognized (Bassi, 2007). One general downside of Python is execution speed. The execution bottleneck of _inmembrane_ lies is in the text processing, where most of the calculations are delegated to the dependent programs that _inmembrane_ forks. As text-processing and system forking are areas in which Python outperforms languages such as Java, there are no performance issues with Python.

Using a modern scripting language results in much cleaner code, where the advantages arise mostly from the use of standard dynamic language features in Python, which otherwise would require the creation of large complicated objects in Java. Another great advantage of scripting languages is portability, where the source code itself is the executable. The source code can be executed directly without any compilation step. Modification of the source-code directly modifies the program. This dramatically simplifies the process of extending _inmembrane_ for other purposes.

### Simple Data Structures 

One unfortunate trend in bioinformatics software is the overuse of object-oriented programming (OOP). In the recommended Enterprise Java style of programming, as used in _SurfG+_, objects are created through several layers of abstract
classes. Each field in an object needs to be minutely specified. To change a field, there are at least 6 places in 3 different files where the code that needs to be changed, which severely restricts the ease of modification. Whilst this level of hierarchy is useful in programs that have highly interacting data-structures, this is not the case for _inmembrane_ and would otherwise add unneeded levels of complexity. 

Python provides an incredibly powerful standard data structure in the form of dictionaries. Dictionaries consists of a set of key-value pairs, where keys and values can be any type of data structure - strings, integers, floats, or even other dictionaries. In _inmembrane_, the program data is represented with a flat dictionary called `protein`. Let's say our FASTA file contains the mouse hemeglobin gene with the ID  `'MOUSE_HEME'`. The properties of `MOUSE_HEME` would then be found in `protein['MOUSE_HEME']`, which is itself a dictionary. `protein['MOUSE_HEME']` contains any arbitary number of different properties, also accessed as key-value pairs. For instance, the sequence length of the `'MOUSE_HEME'` sequence would be stored in `protein['MOUSE_HEME']['sequence_length']`. This data structure can capture the results of any potential bioinformatic analysis, where new properties are added to `protein` on the fly. The use of a dynamic flat dictionary avoids much of the boilerplate involved with an OOP style programming.

### Parsing code is particularly simple

If we use a dictionary to represent our data structure, then the main work in _inmembrane_ of running other programs and processing their text output can be encapsulated into a simple function. 

Using the processing of _SignalP_ as an example, we define a function `run_signalp(protein)` which takes the main protein data structure as input. The function runs SignalP, and then parses the text output. Text processing is very easy to write in Python and the processing for _SignalP_ can be done in about 10 lines of code. As `run_signalp`  cycles through the text output of _SignalP_, then for each protein, say `'MOUSE_HEME'`, if a secretion signal is found, a new property is added: `protein['MOUSE_HEME']['is_signalp'] == True`. 

We can thus abstract the main program loop as running a series of functions of the generic form `run_program(protein)`. This provides a simple API to extend new modules by simply adding new functions that annotates the `protein` dictionary. Depending on the details of how programs are executed, and the format of the output, the functions can get a bit more complicated than `run_signalp`. 

We have exploited other examples of dynamic programming in _inmembrane_ to produce, not only terser code, but clear code. For instance, when using _HMMER_ to match sequence profiles of peptidoglcan binding domains, there is no hard-coding of the sequences in program. Rather, the code dynamically searches the `hmm_profiles` directory for profiles, and iterates the search for each FASTA file across the profiles. This is a particularly robust design, as new profiles can be processed by simply adding them to the directory.

### Web-interface for FASTA file

Hi Andrew!

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

Robel Y. Kahsay1, Guang Gao1 and Li Liao1. An improved hidden Markov model for transmembrane protein detection and topology prediction and its applications to complete genomes (2005) Bioinformatics 21: 1853-1858.

N.Y. Yu, J.R. Wagner, M.R. Laird, G. Melli, S. Rey, R. Lo, P. Dao, S.C. Sahinalp, M. Ester, L.J. Foster, F.S.L. Brinkman (2010) PSORTb 3.0: Improved protein subcellular localization prediction with refined localization subcategories and predictive capabilities for all prokaryotes, Bioinformatics 26(13):1608-1615 <http://dx.doi.org/10.1093/bioinformatics/btq249>

