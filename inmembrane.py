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

import sys
import os
import time
import StringIO
import math
import glob
import re
import shutil
from optparse import OptionParser
import twill
from twill.commands import find, formfile, follow, fv, go, show, \
                             showforms, showlinks, submit

from helpers import *


# Load all modules found in plugins dynamically 
# Each module should have the structure where there is only 
# one main function with the same name of the function
# and takes two parameters:
#   mymodule.mymodule(params, proteins)
for plugin in glob.glob('plugins/*.py'):
  if "__init__" in plugin:
    continue
  plugin_name = os.path.basename(plugin)[:-3]
  exec('from plugins.%s import * ' % plugin_name)


# when True, dumps lots of raw info to stdout to help debugging
__DEBUG__ = False

default_params_str = """{
  'fasta': '',
  'output': '',
  'out_dir': '',
  'organism': 'gram+',
  'signalp4_bin': 'signalp',
  'lipop1_bin': 'LipoP',
  'tmhmm_bin': 'tmhmm',
  'helix_programs': ['tmhmm'],
#' helix_programs': ['tmhmm', 'memsat3'],
  'barrel_programs': ['bomp', 'tmbeta'],
# 'barrel_programs': ['bomp', 'tmbeta', 'tmbhunt'],
  'bomp_cutoff': 1,
  'tmbhunt_cutoff': 0.5,
  'memsat3_bin': 'runmemsat',
  'hmmsearch3_bin': 'hmmsearch',
  'hmm_profiles_dir': '%(hmm_profiles)s',
  'hmm_evalue_max': 0.1,
  'hmm_score_min': 10,
  'terminal_exposed_loop_min': 50,
  'internal_exposed_loop_min': 100,
}
"""


def abs_config_fname():
  module_dir = os.path.abspath(os.path.dirname(__file__))
  return os.path.join(module_dir, 'inmembrane.config')


def get_params():
  config = abs_config_fname()
  if not os.path.isfile(config):
    error_output("# Couldn't find inmembrane.config file")
    error_output("# So, will generate a default config " \
                  + os.path.abspath(config))
    abs_hmm_profiles = os.path.join(module_dir, 'hmm_profiles')
    default_str = default_params_str % \
        { 'hmm_profiles': abs_hmm_profiles }
    open('inmembrane.config', 'w').write(default_str)
  else:
    error_output("# Loading existing inmembrane.config")
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

  fasta = "input.fasta"
  shutil.copy(params['fasta'], os.path.join(base_dir, fasta))
  params['fasta'] = fasta

  shutil.copy(abs_config_fname(), os.path.join(base_dir, config_file))

  os.chdir(base_dir)


def create_protein_data_structure(fasta):
  """
  From a FASTA file, creates a dictionary using the ID of each sequence
  in the fasta file. Also returns a list of ID's in the same order as
  that in the file.
  """
  prot_ids = []
  prot_id = None
  proteins = {}
  for l in open(fasta):
    if l.startswith(">"):
      prot_id, name = parse_fasta_header(l)
      prot_ids.append(prot_id)
      proteins[prot_id] = {
        'seq':"",
        'name':name,
      }
      continue
    if prot_id is not None:
      words = l.split()
      if words:
        proteins[prot_id]['seq'] += words[0]
  return prot_ids, proteins


def chop_nterminal_peptide(protein, i_cut):
  protein['sequence_length'] -= i_cut
  for prop in protein:
    if '_loops' in prop or '_helices' in prop:
      loops = protein[prop]
      for i in range(len(loops)):
        j, k = loops[i]
        loops[i] = (j - i_cut, k - i_cut)
  for prop in protein:
    if '_loops' in prop or '_helices' in prop:
      loops = protein[prop]
      for i in reversed(range(len(loops))):
        j, k = loops[i]
        # tests if this loop has been cut out
        if j<=0 and k<=0:
          del loops[i]
        # otherewise, neg value means loop is at the new N-terminal
        elif j<=0 and k>0:
          loops[i] = (1, k)


def eval_surface_exposed_loop(
    sequence_length, n_transmembrane_region, outer_loops, 
    terminal_exposed_loop_min, internal_exposed_loop_min):
    
  if n_transmembrane_region == 0:
    # treat protein as one entire exposed loop
    return sequence_length >= terminal_exposed_loop_min

  if not outer_loops:
    return False

  loop_len = lambda loop: abs(loop[1]-loop[0]) + 1

  # if the N-terminal loop sticks outside
  if outer_loops[0][0] == 1:
    nterminal_loop = outer_loops[0]
    del outer_loops[0]
    if loop_len(nterminal_loop) >= terminal_exposed_loop_min:
      return True

  # if the C-terminal loop sticks outside
  if outer_loops:
    if outer_loops[-1][-1] == sequence_length:
      cterminal_loop = outer_loops[-1]
      del outer_loops[-1]
      if loop_len(cterminal_loop) >= terminal_exposed_loop_min:
        return True

  # test remaining outer loops for length
  for loop in outer_loops:
    if loop_len(loop) >= internal_exposed_loop_min:
      return True

  return False


def predict_surface_exposure(params, protein):

  def sequence_length(protein):
    return protein['sequence_length']
    
  def has_tm_helix(protein):
    for program in params['helix_programs']:
      if dict_get(protein, '%s_helices' % program):
        return True
    return False

  def has_surface_exposed_loop(protein):
    for program in params['helix_programs']:
      if eval_surface_exposed_loop(
          protein['sequence_length'], 
          len(protein['%s_helices' % (program)]), 
          protein['%s_outer_loops' % (program)], 
          params['terminal_exposed_loop_min'], 
          params['internal_exposed_loop_min']):
        return True
    return False

  terminal_exposed_loop_min = \
      params['terminal_exposed_loop_min']

  is_hmm_profile_match = dict_get(protein, 'hmmsearch')
  is_lipop = dict_get(protein, 'is_lipop')
  if is_lipop:
    i_lipop_cut = protein['lipop_cleave_position']
  is_signalp = dict_get(protein, 'is_signalp')
  if is_signalp:
    i_signalp_cut = protein['signalp_cleave_position']

  details = ""
  if is_hmm_profile_match:
    details += "hmm(%s);" % protein['hmmsearch'][0]
  if is_lipop: 
    details += "lipop;"
  if is_signalp:
    details += "signalp;"
  for program in params['helix_programs']:
    if has_tm_helix(protein):
      n = len(protein['%s_helices' % program])
      details += program + "(%d);" % n

  if is_lipop: 
    chop_nterminal_peptide(protein, i_lipop_cut)
  elif is_signalp:
    chop_nterminal_peptide(protein, i_signalp_cut)

  if is_hmm_profile_match:
    category =  "PSE"
  elif has_tm_helix(protein):
    if has_surface_exposed_loop(protein):
      category = "PSE"
    else:
      category = "MEMBRANE"
  else:
    if is_lipop:
      # whole protein considered outer terminal loop
      if sequence_length(protein) < terminal_exposed_loop_min:
        category = "MEMBRANE"
      else:
        category = "PSE"
    elif is_signalp:
      category = "SECRETED"
    else:
      category = "CYTOPLASM"

  return details, category


def identify_pse_proteins(params):
  prot_ids, proteins = create_protein_data_structure(params['fasta'])

  features = [signalp4, lipop1, hmmsearch3]
  if dict_get(params, 'helix_programs'):
    if 'tmhmm' in params['helix_programs']:
      features.append(tmhmm)
    if 'memsat3' in params['helix_programs']:
      features.append(memsat3)
  if dict_get(params, 'barrel_programs'):
    if 'tmbhunt' in params['barrel_programs']:
      features.append(tmbhunt_web)
    if 'bomp' in params['barrel_programs']:
      features.append(bomp_web)
  for extract_protein_feature in features:
    extract_protein_feature(params, proteins)

  for prot_id in prot_ids:
    details, category = \
        predict_surface_exposure(params, proteins[prot_id])
    if details.endswith(';'):
      details = details[:-1]
    if details is '':
      details = "."
    proteins[prot_id]['details'] = details
    proteins[prot_id]['category'] = category
  
  for prot_id in prot_ids:
    protein = proteins[prot_id]
    print '%-15s ,  %-13s , %-50s , "%s"' % \
        (prot_id, 
         protein['category'], 
         protein['details'],
         protein['name'][:60])

  return prot_ids, proteins


def predict_surface_exposure_barrel(params, protein):
  # TODO: This is a placeholder for a function which will do something
  #       similar to predict_surface_exposure, but focussed on inferring 
  #       outer membrane beta barrel topology.
  #       Essentially, we should:
  #        * Move through the strand list in reverse.
  #        * Strand annotation alternates 'up' strand and 'down' strand
  #        * Loop annotation (starting with the C-terminal residue) alternates
  #          'inside' and 'outside'.
  #        * If everything is sane, we should finish on a down strand. If not,
  #          consider a rule to make an 'N-terminal up strand' become 
  #          an 'inside loop'
  #        * Sanity check on loop lengths ? 'Outside' loops should be on average
  #          longer than non-terminal 'inside' loops.
  #        * For alternative strand predictors (eg transFold, ProfTMB), which
  #          may specifically label inner and outer loops, we should obviously
  #          use those annotations directly.
  pass


def print_summary_table(proteins):
  counts = {}
  counts["BARREL"] = 0
  for seqid in proteins:
    category = proteins[seqid]['category']
    
    # WIP: greedy barrel annotation
    if (dict_get(proteins[seqid], 'tmbhunt_prob') >= params['tmbhunt_cutoff']) or \
       (dict_get(proteins[seqid], 'bomp') >= params['bomp_cutoff']):
       counts["BARREL"] += 1
    
    if category not in counts:
      counts[category] = 0
    else:
      counts[category] += 1
      
  error_output("# Number of proteins in each class:")
  for c in counts:
    error_output("%-15s %i" % (c, counts[c]))


def dump_results(proteins):
  for i,d in proteins.items():
    error_output("# %s - %s" % (i, proteins[i]['name']))
    for x,y in d.items():
      error_output(`x`+": "+`y`)


def identify_omps(params, stringent=False):
  """
  Identifies outer membrane proteins from gram-negative bacteria.
  
  If stringent=True, all predicted outer membrane barrels must also
  have a predicted signal sequence to be categorized as BARREL.
  """
  
  seqids, proteins = create_protein_data_structure(params['fasta'])

  features = [signalp4, lipop1, hmmsearch3]
  if dict_get(params, 'helix_programs'):
    if 'tmhmm' in params['helix_programs']:
      features.append(tmhmm)
    if 'memsat3' in params['helix_programs']:
      features.append(memsat3)
  if dict_get(params, 'barrel_programs'):
    if 'tmbhunt' in params['barrel_programs']:
      features.append(tmbhunt_web)
    if 'bomp' in params['barrel_programs']:
      features.append(bomp_web)
  for extract_protein_feature in features:
    extract_protein_feature(params, proteins)
  
  for seqid, protein in proteins.items():
    # TODO: this is used for setting 'category', however
    #       we may need to make a gram- OM specific version
    #       (eg, run after strand prediction so we can look at
    #            strand topology, detect long extracellular loops etc) 
    details, category = predict_surface_exposure(params, protein)
    proteins[seqid]['category'] = category
    proteins[seqid]['details'] = details
    
    if stringent:
      if dict_get(protein, 'is_signalp') and \
       ( dict_get(protein, 'bomp') or \
         dict_get(protein, 'tmbhunt') ):
       proteins[seqid]['category'] = 'BARREL'
    else:
      if dict_get(protein, 'bomp') or \
         dict_get(protein, 'tmbhunt'):
         proteins[seqid]['category'] = 'BARREL'
    
  # TMBETA-NET knows to only run on predicted barrels
  if 'tmbeta' in params['barrel_programs']:
    tmbeta_net_web(params, proteins, category='BARREL')

  for seqid in proteins:
    details = proteins[seqid]['details']
    if dict_get(proteins[seqid], 'tmbeta_strands'):
      num_strands = len(proteins[seqid]['tmbeta_strands'])
      details += 'tmbeta(%i)' % (num_strands)
    if details.endswith(';'):
      details = details[:-1]
    if details is '':
      details = "."
    proteins[seqid]['details'] = details
    
  print_summary_table(proteins)
  #dump_results(proteins)

  return seqids, proteins
  
  
class Logger(object):
    def __init__(self, log_fname):
        self.terminal = sys.stdout
        self.log = open(log_fname, 'w')

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)  


def process(params):
  if dict_get(params, 'output'):
    sys.stdout = Logger(params['output'])
  init_output_dir(params)
  if params['organism'] == 'gram+':
    seqids, proteins = identify_pse_proteins(params)
  elif params['organism'] == 'gram-':
    seqids, proteins = identify_omps(params, stringent=False)
  else:
    error_output("You must specify 'gram+' or 'gram-' in inmembrane.config\n")
    

description = """
Inmembrane is a proteome annotation pipeline. It takes 
a FASTA file, then carries out sequential analysis of 
each sequence with a bunch of third-party programs, and 
collates the results.

(c) 2011 Bosco Ho and Andrew Perry
"""

if __name__ == "__main__":
  parser = OptionParser()
  (options, args) = parser.parse_args()
  params = get_params()
  if ('fasta' not in params or not params['fasta']) and not args:
    error_output(description)
    parser.print_help()
    sys.exit(1)
  if 'fasta' not in params or not params['fasta']:
    params['fasta'] = args[0]
  process(params)

