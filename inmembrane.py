#!/usr/bin/env python
# Copyright (c) 2012, Bosco Ho <boscoh@gmail.com>
# Copyright (c) 2012, Andrew Perry <ajperry@pansapiens.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the following 
# conditions are met:

# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
# DAMAGE.

__version__ = "0.9"

import sys
import os
import glob
import shutil
from optparse import OptionParser


from helpers import *

# will load all plugins in the plugins/ directory
from plugins import *

description = """
inmembrane %s (https://github.com/boscoh/inmembrane)

Inmembrane is a proteome annotation pipeline. It takes 
a FASTA file, then carries out sequential analysis of 
each sequence with a bunch of third-party programs, and 
collates the results.

(c) 2011 Bosco Ho and Andrew Perry
""" % (__version__)


# figure out absoulte directory for inmembrane scripts
module_dir = os.path.abspath(os.path.dirname(__file__))


default_params_str = """{
  'fasta': '',
  'csv': '',
  'out_dir': '',
  'protocol': 'gram_pos', # 'gram_neg'
  'signalp4_bin': 'signalp',
  'lipop1_bin': 'LipoP',
  'tmhmm_bin': 'tmhmm',
  'helix_programs': ['tmhmm'],
#' helix_programs': ['tmhmm', 'memsat3'],
  'barrel_programs': ['bomp', 'tmbeta'],
# 'barrel_programs': ['bomp', 'tmbeta', 'tmbhunt'],
  'bomp_clearly_cutoff': 3, # >= to this, always classify as an OM(barrel)
  'bomp_maybe_cutoff': 1, # must also have a signal peptide to be OM(barrel)
  'tmbhunt_clearly_cutoff': 0.95,
  'tmbhunt_maybe_cutoff': 0.5,
  'tmbetadisc_rbf_method': 'aadp', # aa, dp, aadp or pssm
  'memsat3_bin': 'runmemsat',
  'hmmsearch3_bin': 'hmmsearch',
  'hmm_profiles_dir': '%(hmm_profiles)s',
  'hmm_evalue_max': 0.1,
  'hmm_score_min': 10,
  'terminal_exposed_loop_min': 50, # unused in gram_neg protocol
  'internal_exposed_loop_min': 100, # try 30 for gram_neg
}
"""


def get_params():
  """
  Gets the params dictionary that hold all the configuration
  information of the program. This is loaded from 'inmembrane.config'
  which should be found in the same place as the main binary.

  If 'inmembrane.config' is not a found, a default 'inmembrane.config'
  is generated from 'default_params_str'. The config file should
  be edited if the binaries are not available on the path, or have
  different names.
  """
  config = os.path.join(module_dir, 'inmembrane.config')
  if not os.path.isfile(config):
    log_stderr("# Couldn't find inmembrane.config file")
    log_stderr(
         "# So, will generate a default config " + config)
    abs_hmm_profiles = os.path.join(module_dir, 'hmm_profiles')
    default_str = default_params_str % \
        { 'hmm_profiles': abs_hmm_profiles }
    open(config, 'w').write(default_str)
  else:
    log_stderr("# Loading existing inmembrane.config")
  params = eval(open(config).read())
  return params


def init_output_dir(params):
  """
  Creates a directory for all output files and makes it the current 
  working directory. copies the input sequences into it as 'input.fasta'.
  """
  if dict_get(params, 'out_dir'):
    base_dir = params['out_dir']
  else:
    base_dir = '.'.join(os.path.splitext(params['fasta'])[:-1])
    params['out_dir'] = base_dir
  if not os.path.isdir(base_dir):
    os.makedirs(base_dir)

  if not dict_get(params, 'csv'):
    basename = '.'.join(os.path.splitext(params['fasta'])[:-1])
    params['csv'] = basename + '.csv'
  params['csv'] = os.path.abspath(params['csv'])

  params['citations'] = os.path.join(params['out_dir'], 'citations.txt')
  params['citations'] = os.path.abspath(params['citations'])
  
  fasta = "input.fasta"
  shutil.copy(params['fasta'], os.path.join(base_dir, fasta))
  params['fasta'] = fasta

  shutil.copy(os.path.join(module_dir, 'inmembrane.config'), base_dir)

  os.chdir(base_dir)


def import_protocol_python(params):
  """
  Some python magic that loads the desired protocol file
  encoded in the string 'params['protocol'] as a python file
  with the internal variable name 'protocol'. An appropriate
  python command is generated that is to be processed by
  the 'exec' function.
  """
  protocol_py = os.path.join(module_dir, 'protocols', params['protocol']+'.py')
  if not os.path.isfile(protocol_py):
    raise IOError("Couldn't find protcols/" + protocol_py)
  return 'import protocols.%s as protocol' % (params['protocol'])


def process(params):
  """
  Main program loop. Triggers the 'protocol' found in the params
  to annotate all proteins give the list of annotations needed by
  'protocol'. Then outputs to screen and a .csv file.
  """
  # initializations
  exec(import_protocol_python(params))
  init_output_dir(params)
  seqids, proteins = create_proteins_dict(params['fasta'])

  # TODO: ideally this loop needs to be run within the protocol,
  #       since for some protocols not all plugins
  #       will be run for every sequence, conditional
  #       on the outcome of a previous analysis
  #       eg. protocol.run(params, proteins)

  # annotates with external binaries as found in plugins
  for plugin_str in protocol.get_annotations(params):
    plugin = eval(plugin_str)
    plugin.annotate(params, proteins)
  
  # do protocol analysis on the results of the annotations
  for seqid in seqids:
    protein = proteins[seqid]
    protocol.post_process_protein(params, protein)
    log_stdout(protocol.protein_output_line(seqid, proteins))

  # print a summary table of classifications to stderr
  log_stderr(protocol.summary_table(params, proteins))

  # always write to biologist-friendly csv file
  f = open(params['csv'], 'w')
  for seqid in seqids:
    f.write(protocol.protein_csv_line(seqid, proteins))
  f.close()
  log_stderr("\n")
  log_stderr("Output written to %s" % (params['csv']))

  # TODO: citations for specific HMMs (PFAM etc ?)
  
  # write citations to a file and gracefully deal with plugins
  # without a citation defined
  import codecs
  import textwrap
  f = codecs.open(params['citations'], mode='w', encoding='utf-8')
  programs_used = []
  for program in protocol.get_annotations(params):
    plugin = eval(program)
    try:
      f.write(plugin.citation['name']+":\n")
      f.write(textwrap.fill(plugin.citation['ref']))
    except AttributeError:
      f.write("%s - no citation provided." % program)
    f.write("\n\n")
    try:
      programs_used.append(plugin.citation['name'])
    except AttributeError, KeyError:
      programs_used.append(program)
    
  f.close()
  log_stderr("\n")
  log_stderr("This run used %s." % (", ".join(programs_used)) )
  log_stderr("References have been written to %s \n"
             "# - please cite as appropriate." % 
             (params['citations']) )
     
            
if __name__ == "__main__":
  parser = OptionParser()
  (options, args) = parser.parse_args()
  params = get_params()
  if ('fasta' not in params or not params['fasta']) and not args:
    print description
    parser.print_help()
    sys.exit(1)
  if 'fasta' not in params or not params['fasta']:
    params['fasta'] = args[0]
  process(params)

