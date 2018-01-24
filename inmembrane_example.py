#!/usr/bin/env python

params = {
    'fasta': 'examples/surfg/AE004092.fasta',
    'out_dir': '',
    'csv': 'examples/surfg/AE004092.csv',
    'protocol': 'gram_pos',
    'signalp4_bin': 'signalp_scrape_web',
    'signalp4_organism': 'gram+',
    'lipop1_bin': 'lipop_scrape_web',
    'tmhmm_bin': 'tmhmm_scrape_web',
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
