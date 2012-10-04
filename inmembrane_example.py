#!/usr/bin/env python

params = {
  'fasta': 'examples/surfg/AE004092.fasta',
  'out_dir': '',
  'csv': 'examples/surfg/AE004092.csv',
  'protocol': 'gram_pos',
  'signalp4_bin': 'signalp',
  'signalp4_organism': 'gram+',
  'lipop1_bin': 'LipoP',
  'tmhmm_bin': 'tmhmm',
  'helix_programs': ['tmhmm'],
  'memsat3_bin': 'runmemsat',
  'hmmsearch3_bin': 'hmmsearch',
  'hmm_profiles_dir': 'protocols/gram_pos_profiles',
  'hmm_evalue_max': 0.1,
  'hmm_score_min': 10,
  'terminal_exposed_loop_min': 50,
  'internal_exposed_loop_min': 100,
}

import inmembrane
inmembrane.process(params)
