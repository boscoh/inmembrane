

# _inmembrane_, a bioinformatic workflow for annotation of bacterial cell-surface proteomes

Here's the situation: you have a bunch of protein sequences, probably in a fasta form, and there's not information about them. So you think that throwing them at a bunch of sequence prediction programs to "guess" what they might do, or how they might act in your organism.

## Combining analyses into a workflow

A more sophisticated analysis combines analyses from different sequence analysis program. In the old days, you would have to download binaries to run the analysis for each analysis. As anyone who's had experience with such binaries know, it's a difficult process, as binaries are provided for different platforms, and there is a whole bunch of different installation procedures involved with all sorts of different programming models.

Fortunately, tody, many sequence-based predictors provide excellent web-interfaces - just do it on the web.However, if you want to do a more involved analysis, then you have to manually filter results from one website to another. 

## _inmembrane_ is a workflow

Here's an example of a workflow - combining predictions of lipo-protein signals, secretion signals, and transmembrane helices - we can identify proteins that are embedded in the membranes, and then subsequently analyze whether they are large enough to poke out of the cell-membrane in the cell.

This requires a combination of several type of prediction programs. 


(c) 2012, Andrew J. Perry & Bosco K. Ho



