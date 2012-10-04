# -*- coding: utf-8 -*-
from inmembrane.helpers import run, parse_fasta_header

citation = {'ref': u"Petersen TN, Brunak S, von Heijne G, Nielsen H. "
                   u"SignalP 4.0: discriminating signal peptides from "
                   u"transmembrane regions. Nature methods 2011 "
                   u"Jan;8(10):785-6. \n"
                   u"<http://dx.doi.org/10.1038/nmeth.1701>",
            'name': "SignalP 4.0"
           }

def annotate(params, proteins):
  for seqid in proteins:
    proteins[seqid]['is_signalp'] = False
    proteins[seqid]['signalp_cleave_position'] = None

  signalp4_out = 'signalp.out'
  cmd = '%(signalp4_bin)s -t %(signalp4_organism)s  %(fasta)s' % \
             params
  run(cmd, signalp4_out)

  for line in open(signalp4_out):
    if line.startswith("#"):
      continue
    words = line.split()
    seqid = parse_fasta_header(words[0])[0]
    proteins[seqid]['signalp_cleave_position'] = int(words[4])
    if (words[9] == "Y"):
      proteins[seqid]['is_signalp'] = True

  return proteins
    
