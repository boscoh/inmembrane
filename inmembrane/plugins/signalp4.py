# -*- coding: utf-8 -*-
from inmembrane.helpers import run, parse_fasta_header

citation = {'ref': u"Petersen TN, Brunak S, von Heijne G, Nielsen H. "
                   u"SignalP 4.0: discriminating signal peptides from "
                   u"transmembrane regions. Nature methods 2011 "
                   u"Jan;8(10):785-6. \n"
                   u"<http://dx.doi.org/10.1038/nmeth.1701>",
            'name': "SignalP 4.0"
            }


def parse_signalp(signalp4_lines, proteins, id_mapping=None):
    if id_mapping is None:
        id_mapping = []

    past_preamble = False
    for line in signalp4_lines:
        if line.startswith("#"):
            past_preamble = True
            continue
        if not past_preamble and line.strip() == '':  # skip empty lines
            continue
        if past_preamble:
            if line.strip() == '':
                # in the case of web output of concatenated signalp output
                # files, an empty line after preamble means we have finished all
                # 'result' lines for that section
                past_preamble = False
                continue
            words = line.split()
            seqid = parse_fasta_header(words[0])[0]
            if id_mapping:
                seqid = id_mapping[seqid]
            proteins[seqid]['signalp_cleave_position'] = int(words[4])
            proteins[seqid]['is_signalp'] = (words[9] == 'Y')

    return proteins


def annotate(params, proteins):
    for seqid in proteins:
        proteins[seqid]['is_signalp'] = False
        proteins[seqid]['signalp_cleave_position'] = None

    signalp4_out = 'signalp.out'
    cmd = '%(signalp4_bin)s -t %(signalp4_organism)s  %(fasta)s' % \
          params
    run(cmd, signalp4_out)

    with open(signalp4_out) as signalp4_text:
        proteins = parse_signalp(signalp4_text, proteins)

    return proteins
