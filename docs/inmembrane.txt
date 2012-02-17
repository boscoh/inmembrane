# Annotation of proteomes: bioinformatic workflow with __inmembrane__, a surface membrane predictor for gram+ bacteria

## Abstract 

__inmembrane__ is a tool to predict surface membrane proteins from a set of sequences. It is  intended to be a direct replacement for surfg+, which designed the protocol for surface membrane analysis. As well __inmembrane__ demonstrates the dramatic simplification in code through using a modern scripting language, which is particularly suited to such analysis of such bioinformatic problems. This leads to a program that is much more robust and amenable to modification and extensibility. This provides an example of writing such glue programs for bioinformatic analysis. 

## Introduction

A common task in bioinformatics is to integrate the results of protein prediction programs to deduce complex properties of proteins. As sequencing techniques improve, the high volume of data will require a greater capacity to process such sequences in an automated way to provide high-level analysis.

One example is the program SurfG+, which is designed to predict proteins that are exposed on the surface of Gram+ bacteria. SurfG+ is a JAVA program that carries out batch processing of several standard bioinformatic tools to specifically identify Gram+ bacterial proteins that may be exposed out of the peptidoglycan layer of the bacterium.These predictions are meant to identify a set of proteins that would be amenable to cell shaving experiments. SurfG+ itself does not carry out any extensive analysis, but rather relies on a transmembrane helix predictor (THMMOD), a secretion signal predictor (SignalP), a lipoprotein signal predictior (LipoP) and a sequence alignment for protein profiles (HMMER). All these programs function as  Linux command tools that take FASTA sequences as input and generates output as formatted text. 

Nevertheless, __SurfG+__ suffers several problems that plague much bioinformatic software. The most egregrious being that the program is not generally available anymore. We were able to find a poorly-managed source-code repository but the we were not able to get the program to work. With the source code, it was possible to deduce the basic algorithm of the program, however, the software architecture of SurfG+ mdade it impossible to debug properly. Although it may appear at first, to be a redundant effort, there are several compelling reasons for rewriting SurfG+. First, the URL mentioned in the original reference no longer exists, even though the paper was only published in 2009. Second, after a google search, we managed to find a compiled version of Surfg+ but could not get it to work under Mac OSX or Ubuntu Linux. Third, after more searching, we managed to find a copy of the source-code. Unfortunately, due to the architecture of the program, it would have been prohibitively time-consuming to debug the program to work on our systems.  

However, the core of __SurfG+__ is rather simple and consequently, we decided to write a tool to reproduce the functionality of __SurfG+__. In the process of improving the architecture of SurfG+ in writing __inmembrane__, we have generated a better template for this class of bioinformatic software. One example is that from a JAVA source of 700K, we have reduced it down 8K of PYTHON code, where the terseness of the code greatly facilitates modifiability. Here, we discuss the issues involved in writing robust bioinformatic source code.

## Discussion


### Public Open Source Repository

Perhaps the single most important step is to distribute our code on an open-source repository, Github. Our code is public, and is presented in best-of-breed. There is excellent code browsing, with history. As well, convenient downloads are automatically generated, with well-defined URL links. These are generally not provided in typical academic releases of software including the dreaded login and register with some lab manager. As well, Github provides download statistics.

More importantly, the code will be guaranteed to stay in existence, in a way that most labs cannot provide for. As well, Github provides excellent facilities to fork a project. That is, if you were to come across a project that did something useful but had been abandoned a while ago, it is almost trivial to create a duplicate and that make changes. This is the best guarantee for orphan code

### Program Workflow 

The heart of SurfG+ is actually quite simple. SurfG+ is a wrapper that takes FASTA sequences, and runs a number of bioinformatic programs (HMMER, LipoP, SignalP and THMMOD) that also take FASTA sequences as input, parses the result of these programs, and then generates a text report. The bulk of the computation occurs in the other programs, where most of the work is in parsing the text output. A very small amount of analysis is done at the end to produce the text output.

### Scripting Languages 

However, one might argue that one big advantage of Java over scripting languages in speed, but there is very little computation done in the workflow. Indeed, the bulk of the program is spent parsing text, and running other programs. These are two areas where scripting languages such as Python outperform languages such as Java. Consequently, writing __inmembrane__ in a modern scripting language results in a vastly simpler program that is easier to modify, where __inmembrane__ is smaller by almost two orders of magnitue in code size. The simplification arises mostly from the use of standard language features that have to be written from scratch in Java. 

The other major saving is the use of powerful and flexible native data structures in Python. Due to the dynamic nature of Python data structures, very little code is needed to add extra modules to the pipeline as we will demonstrate below.

Another great advantage of scripting languages is that the source code is the executable. As long as Python is installed on the system, the actual script itself is all that is needed to run the program. This makes it very easy to deploy on your system. Modifications to the program is instantly applied. However, there is one last caveat as regards the cross-platform of __inmembrane__. As some of the executables it depends on only have Linux executables, it will only be fully operational on a Linux platform. 

Something about intermediate output files. 


### Data Structures

One of the limitations of an object-oriented style of programming in a static-typed system such as Java is that data structures must be minutely specified. Following the recommended Enterprise Java style of programming, abstract objects are created, resulting in a massive data object with a lot of methods that are used only once, if at all. In comparison, in SurfG+, every separate progam requires the creation of data object, each with a getter and setter for every variable, as well as matching structures in the file reading and analysis. For every field, there are at least 6 places in 3 different files where the code that needs to be changed. 

For a small program such as SurfG+, the abstraction results in a dispersal of code that makes it excceedingly difficult to read and to modify. One of the liberating features of modern scripting languages is the provision of powerful data structures in the form of dynamic hash tables, or dictionaries in python. A dictionary consists of a set of key-value pairs, where keys and values can be any type of data structure. 

The basic data structure in __inmembrane__ is a single dictionary called `protein` consisting of key-value pairs where the key is the FASTA ID in the FASTA file, and the value is a dictionary of properties. Given a protein 'ECOLG_HEME', the properties of the protein are `protein['ECOLG_HEME']`. The properties, say the `'sequence_length'`, can simply be accessed by `protein['ECOLG_HEME']['sequence_length']`.

Since Python allows dynamic data modification, each module  simply add properties to the dictionary of properties, as needeed. There is no need to design a dedicated data object where the data structure is prespecified.


### Text-processing function as annotation of dictionary

Whilst in a more complex program, allowing dyanmic modification of the data structure could lead to obfuscated code, in this particular example, we generate clearer code. The reason is that it is very easy to write parsing code in Python. The resultant code is so short, we can embed the running of a program, say SignalP, the parsing of the output text file, and the modification of the main data structure in a very short function.

Typically, each module is embodied in a function, for example `run_signalp(protein)` which takes the main protein data structure, runs SignalP, reads the output, extracts from the output whether a secretion signal is found for all the proteins in `protein`, and annotates the `protein` structure, by assigning new properties to each protein `protein['ECOLG_HEME']['is_signalp'] == True`. 

In Python, the built in grammar is very powerful, and text processing can be expressed in very terse, yet readable. Indeed, you cannot abstract the text-processing any better - trying to create an object that is flexible enough to handle different text-output actually results in more complicated code. It is a fact of life that new versions of these programs will result in some incompatible text structure so that the formatting will break. We take the philosophy that it is better to expose the text parsing, so that with newer versions, it is easier to just rewrite the text parsing, rather than try to anticipate any future changes by writing complex Parsing objects.

Taking this functional approach - where we encapsulates the analysis of each program in a simple form `run_program(protein)`, we can thus abstract the pipeline, as a series of functions with exactly the same interface. This defines a very simple API to extend new modules. Simply write a new function that annotates the data structure.

As such, we have shown to extend the workflow of __inmembrane__ with respect to SurfG by adding an alternative to __THMHMM__, __MEMSAT3__, which produces a much more complex operation, in that __MEMSAT3__ can only take FASTA files with only one sequence. As well, we have added a beta-barrel analysis, which uses a web interface. We have this to show how one might carry out these kinds of extensions.


### Other metaprogramming curios

There are also opportunities to write more flexible code. For instance in the search of HMMER profiles to identify sequence motifs of proteins that directly bind to the peptidoglycan wall, in __inmembrane__, the program searches directly in the directory for profiles, and executes the HMM search as it finds it. In contrast, in SurfG+, each individual profile requires the creation of a separate object leading to a proliferation of classes for each individual sequence motif. The main loop requires a massive switch statement triggering each object.

In python, config files can be coded as a simple dictionary object that is read directly into the program. The `eval` function in Python is used to convert the config file directly into a live `params` dictionary in Python. This allows the use of Python constructs in the config file.


## Test set

unit tests
test with data


## Conclusion

The clarity of Python code better reflects the work-flow described in Surfg+. Python code is often considered the nearest code to pseudo-code and as such, we can use actual code to illustrate the work-flow. 

## References

Barinov A, Loux V, Hammani A, Nicolas P, Langella P, et al. (2009) Prediction of surface exposed proteins in Streptococcus pyogenes, with a potential application to other Gram-positive bacteria. __Proteomics__ 9: 61-73. <http://dx.doi.org/10.1002/pmic.200800195>


