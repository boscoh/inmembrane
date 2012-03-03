# __inmembrane__, a bioinformatic workflow for annotation of bacterial cell surface proteomes

Andrew Perry and Bosco Ho  
_Department of Biochemistry, Monash University, Melbourne, Australia_

## Abstract 

_inmembrane_ is a tool to predict the surface exposed regions of membrane proteins in sets of bacterial protein sequences. It is intended to be a direct replacement for SurfG+, which implemented a protocol for predicting extracellular regions of polypeptide in Gram positive bacterial proteomes. _inmembrane_ provides a more accessible and transparent code base through use of a modern scripting language. This leads to a program that is much more amenable to modification and extensibility, and provides an example of writing such glue programs for bioinformatic analysis. The program is hosted on a public open-source repository http://github.com/boscoh/inmembrane.

## Introduction

A common task in bioinformatics is to integrate the results of protein prediction programs to deduce complex properties of proteins. As sequencing techniques improve, the burgeoning volume of data will require a greater capacity to analyze such sequences in an automated way.

In studies of membrane proteomes, quick annotation of an experimentally detected set of the proteins can help detect sequences of unexpected localization, and can alert researchers to possible contamination from other subcellular fractions. Ultimately, a concise summary of the properties of the detect membrane proteins in a particular proteomic set allows meaningful comparisons between different bacterial strains, species, and their responses in membrane remodelling to host and enviromental challenges.

One example of this type of annotation is the program SurfG+, which is designed to predict proteins that are exposed on the surface of Gram+ bacteria. SurfG+ is a Java program that carries out batch processing of several standard bioinformatic tools to specifically identify Gram+ bacterial proteins that may be exposed out of the peptidoglycan layer of the bacterium. These predictions are intended to identify a set of proteins that would be amenable to cell surface protease shaving experiments. SurfG+ itself does not carry out any extensive analysis, but rather relies on a transmembrane helix predictor (_TMMOD_), a secretion signal predictor (_SignalP_), a lipoprotein signal predictior (_LipoP_) and a sequence alignment for protein profiles (_HMMER_). All these programs function as Linux command tools that take FASTA sequences as input and generates output as formatted text.

Nevertheless, _SurfG+_ suffers several problems that plague much bioinformatic software. The most egregrious being that the program is not generally available anymore. The URL mentioned in the original reference no longer exists, even though the paper was only published in 2009. We were able to find a poorly-managed source-code repository but we could not get the program to work. The program was written in the enterprise Java style of programming, making it extremely difficult to read, and almost impossible to debug. 

Since the core algorithm in _SurfG+_ is relatively straightforward, we decided to write _inmembrane_ to replicate the functionality of _SurfG+_, but in a modern scripting language. This lead to a major simplifiction of the code base - from the original Java source of 700K, we ended it with 10K of Python code. Although none of the techniques used in _inmembrane_ is particularly new, they may be unfamiliar to the bioinformatic community. Here, we present the ideas that lead to, not only a smaller but a clearer code case, which substantially makes it easier for others to reuse and repurpose the code in _inmembrane_.

## Discussion

### Public Open Source Repository

Perhaps the single most important step is to distribute our code on an open-source repository, Github. There is excellent code browsing, with history. As well, convenient downloads are automatically generated, with well-defined URL links. These are generally not provided in typical academic releases of software including the dreaded login and registeration pages. Github provides download statistics, which are useful as a metric to determine the popularity and broader impact of the software.

More importantly, the source code has a higher probability of remaining accessible in the long term, something that historically most academic labs have shown they cannot provide through in-house hosting. As well, Github provides excellent facilities to fork a project. That is, if you were to come across a project that did something useful but had been abandoned a while ago, it is almost trivial to create a duplicate and that make changes. This is the best guarantee for orphan code.

The algorithms used here are not difficult, and there is no reason to hide the program behind a hybrid academic/commercial license. We have licensed under the BSD license, which liberally allows reuse of the code in any shape or form.

### Program Workflow 

The heart of _inmembrane_ is quite simple. It is a wrapper that takes FASTA sequences and sequentially runs a number of sequence analysis programs (HMMER, LipoP, SignalP and TMMOD) that share the common feature of taking FASTA sequences as input, and generating text output. The bulk of the computation in _inmembrane_ lies in the parsing of the result of these programs. A small amount of analysis is done at the end to integrate the produce the final text output. 

As _inmembrane_ integrates the intermediate output of a large number of program, there is a lot of opportunity for failure. _inmembrane_ saves all intermediate output in a results folder, which allows expert users to diagnose potential problems with the dependencies.

As some of the dependent executables of _inmembrane_ are only available on Linux, this unfortunately restricts _inmembrane_ to be fully operational only on the Linux platform. 


### Scripting Languages 

The virtues of Python as a language for solving problems in life science research have been previously recognized (Bassi, 2007). One general downside of Python is execution speed. However, the bottleneck in _inmembrane_ is in text parsing, most of the calculations are delegated to the dependent programs that _inmembrane_ forks. As text-processing and system forking are areas in which Python outperforms other languages such as Java, there are no performance issues in the use of Python.

Using a modern scripting language results in much cleaner code, where the advantages arises mostly from the use of standard dynamic language features in Python, which otherwise would require the creation of large complicated objects in Java. Another great advantage of scripting languages is portability, where the source code itself is the executable. The source code can be executed directly without any compilation step. Modification of the source-code directly modifies the program. This dramatically simplifies the process of extending _inmembrane_ for other purposes.

### Simple Data Structures 

One unfortunate trend in bioinformatics software is the overuse of object-oriented programming (OOP). In the recommended Enterprise Java style of programming, as used in _SurfG+_, objects are created through several layers of abstract layers. Each field in an object needs to be minutely specified. To change a field, there are at least 6 places in 3 different files where the code that needs to be changed, which severly restricts the ease of modification.

Whilst this level of hierarchy is useful in programs that have highly interacting data-structures, there is hardly any interaction between the different fields in _inmembrane_ and programming multiple levels of object hiearchy just adds unneeded levels of complexity. 

In _inmembrane_, the program data is represented with a standard dictionary called `protein`. Dictionaries consists of a set of key-value pairs, where keys and values can be any type of data structure. For instance, the FASTA ID of the mouse hemeglobin gene might be `'MOUSE_HEME'`. The properties of `MOUSE_HEME` is thus found in `protein['MOUSE_HEME']`. 

`protein['MOUSE_HEME']` is itself a dictionary, which contain any arbitary number of different properties, also accessed as key-value pairs. For instance, to get the sequence length of the `'MOUSE_HEME'` sequence, you type `protein['MOUSE_HEME']['sequence_length']`
     
New properties are added to `protein` on the fly, and thus avoids the boilerplate involved with an OOP style programming.

### Parsing code is particularly simple

If we use a dictionary to represent our data structure, then the main work in _inmembrane_ of running other programs and processing the text output can be encapsulated into a simple function. 

Using the processing of _SignalP_ as an example, we define a function `run_signalp(protein)` which takes the main protein data structure as input. The function runs SignalP, and then parses the text output. Text processing is very easy to write in Python and the processing for _SignalP_ can be done in about 10 lines of code. As `run_signalp`  cycles through the text output of _SignalP_, then for each protein, say `'MOUSE_HEME'`, if a secretion signal is found, a new property is added: `protein['MOUSE_HEME']['is_signalp'] == True`. 

We can thus abstract the main program loop as running a series of functions of the generic form `run_program(protein)`. This provides a simple API to extend new modules by simply adding new functions that annotates the `protein` dictionary. Depending on the details of how programs are executed, and the format of the output, the functions can get a bit more complicated than `run_signalp`. 

### Other dynamic programming techniques

We have exploited other examples of dynamic programming in _inmembrane_ to produce, not only terser code, but clear code. For instance, when using _HMMER_ to match sequence profiles of peptidoglcan binding domains, there is no hard-coding of the sequences in program. Rather, the code dynamically searches the `hmm_profiles` directory for profiles, and iterates the search for each FASTA file across the profiles. This is a particularly robust design, as new profiles can be processed by simply adding them to the directory.

In python, config files can be coded as a simple dictionary object that is read directly into the program. The `eval` function in Python is used to convert the config file directly into a live `params` dictionary in Python. This allows the use of Python constructs in the config file.


## Test set

unit tests
test with data


## Conclusion

The clarity of Python code better reflects the work-flow described in Surfg+. Python code is often considered the nearest code to pseudo-code and as such, we can use actual code to illustrate the work-flow. 

## References

Barinov A, Loux V, Hammani A, Nicolas P, Langella P, et al. (2009) Prediction of surface exposed proteins in Streptococcus pyogenes, with a potential application to other Gram-positive bacteria. __Proteomics__ 9: 61-73. <http://dx.doi.org/10.1002/pmic.200800195>

Bassi S (2007) A Primer on Python for Life Science Researchers. __PLoS Comput Biol__ 3(11): e199. <http://dx.doi.org/10.1371/journal.pcbi.0030199>

Jannick Dyrløv Bendtsen, Henrik Nielsen, Gunnar von Heijne and Søren Brunak. (2004) Improved Prediction of Signal Peptides: SignalP 3.0. __J. Mol. Biol.__ 340:783–795.

Robert D. Finn, Jody Clements and Sean R. Eddy. (2011) HMMER web server: interactive sequence similarity searching. __Nucleic Acids Research__ 39:W29–W37.

Agnieszka S. Juncker, Hanni Willenbrock, Gunnar Von Heijne, Søren Brunak, Henrik Nielsen, And Anders Krogh. (2003) Prediction of lipoprotein signal peptides in Gram-negative bacteria. __Protein Science__ 12:1652–1662.

Anders Krogh, Björn Larsson, Gunnar von Heijne and Erik L. L. Sonnhammer (2001) Predicting Transmembrane Protein Topology with a Hidden Markov Model: Application to Complete Genomes. __J. Mol. Biol.__ 305:567-580.

TMMOD


