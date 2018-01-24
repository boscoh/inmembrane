# -*- coding: utf-8 -*-
import os
import re
from inmembrane.helpers import seqid_to_filename, run, write_proteins_fasta

citation = {'ref': u"﻿Jones DT. Improving the accuracy of transmembrane "
                   u"protein topology prediction using evolutionary "
                   u"information. Bioinformatics (Oxford, England) "
                   u"2007 Mar;23(5):538-44. \n"
                   u"<http://dx.doi.org/10.1093/bioinformatics/btl677",
            'name': "MEMSAT3"
            }


def parse_memsat(protein, memsat_out):
    """
    Parses a MEMSAT3 output file for a FASTA sequence and
    annotates the protein dictionary. This is an auxillary function
    called by annotate, and the fields are described there.
    """

    f = open(memsat_out)
    l = f.readline()
    while l:
        l = f.readline()

        if l == "FINAL PREDICTION\n":
            f.readline()
            l = f.readline()
            s = l.split(":")

            while re.match("\d", l[0]):
                tokens = s[1].strip().split()
                tok_offset = 0
                if len(tokens) > 2:
                    tok_offset = 1
                    side_of_membrane_nterminus = tokens[0][
                                                 1:-1]  # 'in' or 'out'
                i = int(tokens[tok_offset].split('-')[0])
                j = int(tokens[tok_offset].split('-')[1])
                protein['memsat3_helices'].append((i, j))
                score = float(tokens[1 + tok_offset][1:-1])
                protein['memsat3_scores'].append(score)
                l = f.readline()
                s = l.split(":")
            f.readline()

            # record inner and outer loops
            inner_loops = protein['memsat3_inner_loops']
            outer_loops = protein['memsat3_outer_loops']
            sequence_length = protein['sequence_length']
            if side_of_membrane_nterminus == 'out':
                loops = outer_loops
            elif side_of_membrane_nterminus == 'in':
                loops = inner_loops
            loop_start = 1

            # figure out helices
            for tm in protein['memsat3_helices']:
                loop_end = tm[0] - 1
                loop = (loop_start, loop_end)
                loops.append((loop_start, loop_end))
                if loops == outer_loops:
                    loops = inner_loops
                else:
                    loops = outer_loops
                loop_start = tm[1] + 1

                # capture C-terminal loop segment
            loop_start = tm[1] + 1
            loop_end = sequence_length
            loops.append((loop_start, loop_end))

    f.close()


def has_transmembrane_in_globmem(globmem_out):
    for l in open(globmem_out):
        if "Your protein is probably not a transmembrane protein" in l:
            return False
    return True


def annotate(params, proteins):
    """
    Runs MEMSAT3 and parses the output files. Takes a standard 'inmembrane'
    params dictionary and a global proteins dictionary which it populates with
    results.

    In the current implementation, this function extracts and feeds sequences to MEMSAT3
    one by one via a temporary file.

    These keys are added to the proteins dictionary:
      - 'memsat3_helices', a list of tuples describing the first and last residue
        number of each transmembrane segment;

      - 'memsat3_scores', a list of confidence scores (floats) for each predicted
        tm segment;

      - 'memsat3_inner_loops', a list of tuples describing the first and last residue
         number of each predicted internal loop segment;

      - 'memsat3_outer_loops', a list of tuples describing the first and last residue
         number of each predicted outer loop segment;
    """

    for seqid in proteins:
        protein = proteins[seqid]

        # initialize the protein data structure
        protein.update({
            'memsat3_scores': [],
            'memsat3_helices': [],
            'memsat3_inner_loops': [],
            'memsat3_outer_loops': []
        })

        # write seq to single fasta file
        single_fasta = seqid_to_filename(seqid) + '.fasta'
        if not os.path.isfile(single_fasta):
            write_proteins_fasta(single_fasta, proteins, [seqid])

        memsat_out = single_fasta.replace('fasta', 'memsat')
        run('%s %s' % (params['memsat3_bin'], single_fasta), memsat_out)

        globmem_out = single_fasta.replace('fasta', 'globmem')
        if has_transmembrane_in_globmem(globmem_out):
            parse_memsat(protein, memsat_out)
