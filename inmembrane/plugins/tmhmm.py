# -*- coding: utf-8 -*-

citation = {'ref': u"Anders Krogh, Björn Larsson, Gunnar von Heijne and Erik "
                   u"L. L. Sonnhammer (2001) Predicting Transmembrane Protein "
                   u"Topology with a Hidden Markov Model: Application to "
                   u"Complete Genomes. J. Mol. Biol. 305:567-580. \n"
                   u"<http://dx.doi.org/10.1006/jmbi.2000.4315>",
            'name': "TMHMM 2.0"
            }

from inmembrane.helpers import run, parse_fasta_header


def annotate(params, proteins):
    """
    Runs THMHMM and parses the output files. Takes a standard 'inmembrane'
    params dictionary and a global proteins dictionary which it populates with
    results.

    In the current implementation, this function extracts and feeds sequences to tmhmm
    one by one via a temporary file.

    These keys are added to the proteins dictionary:
      - 'tmhmm_helices', a list of tuples describing the first and last residue
        number of each transmembrane segment;

      - 'tmhmm_scores', a list of confidence scores (floats) for each predicted
        tm segment;

      - 'tmhmm_inner_loops', a list of tuples describing the first and last residue
         number of each predicted internal loop segment;

      - 'tmhmm_outer_loops', a list of tuples describing the first and last residue
         number of each predicted outer loop segment;
    """

    tmhmm_out = 'tmhmm.out'
    run('%(tmhmm_bin)s %(fasta)s' % params, tmhmm_out)
    return parse_tmhmm(open('tmhmm.out').read(), proteins)


def parse_tmhmm(text, proteins, id_mapping=[]):
    seqid = None
    for i_line, l in enumerate(text.split('\n')):
        if i_line == 0:
            continue
        words = l.split()
        if not words:
            continue

        if l.startswith("#"):
            seqid = parse_fasta_header(words[1])[0]
        else:
            seqid = parse_fasta_header(words[0])[0]
        if seqid is None:
            continue

        # re-map from a safe_seqid to the original seqid
        if id_mapping:
            seqid = id_mapping[seqid]

        # initialize fields in proteins[seqid]
        if 'tmhmm_helices' not in proteins[seqid]:
            proteins[seqid].update({
                'tmhmm_helices': [],
                'tmhmm_inner_loops': [],
                'tmhmm_outer_loops': []
            })

        if 'inside' in l:
            proteins[seqid]['tmhmm_inner_loops'].append(
                (int(words[-2]), int(words[-1])))
        if 'outside' in l:
            proteins[seqid]['tmhmm_outer_loops'].append(
                (int(words[-2]), int(words[-1])))
        if 'TMhelix' in l:
            proteins[seqid]['tmhmm_helices'].append(
                (int(words[-2]), int(words[-1])))

    return proteins
