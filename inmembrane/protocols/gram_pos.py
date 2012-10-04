import os
from inmembrane.helpers import dict_get, eval_surface_exposed_loop, \
                    chop_nterminal_peptide

def get_annotations(params):
  """
  Creates a list of annotation functions required
  by this gram_pos protocol. The main program
  will run the annotation functions of this list,
  mapping the correct functions to the strings.

  As well, the function does some bookeeping on
  params to make sure the 'hmm_profiles_dir' is
  pointing in the right place.
  """
  annotations = []
  
  params['signalp4_organism'] = 'gram+'
  
  if not params['signalp4_bin'] or params['signalp4_bin'] == 'signalp_web':
    annotations += ['signalp_web']
  else:
    annotations += ['signalp4']
  
  if not params['lipop1_bin'] or params['lipop1_bin'] == 'lipop_web':
    annotations += ['lipop_web']
  else:
    annotations += ['lipop1']
    
  annotations += ['hmmsearch3']

  if dict_get(params, 'helix_programs'):
    if 'tmhmm' in params['helix_programs']:
      annotations.append('tmhmm')
    if 'memsat3' in params['helix_programs']:
      annotations.append('memsat3')

  params['hmm_profiles_dir'] = os.path.join(
      os.path.dirname(__file__), 'gram_pos_profiles')

  return annotations


def post_process_protein(params, protein):
  """
  This is the main analysis of the protein, where theprotein
  dictionary should contain all the necessary information
  from the annotations. Thus post_process_protein contain
  can determine the final analysis.
  """

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

  details = []
  if is_hmm_profile_match:
    details += ["hmm(%s)" % "|".join(protein['hmmsearch'])]
  if is_lipop: 
    details += ["lipop"]
  if is_signalp:
    details += ["signalp"]
  for program in params['helix_programs']:
    if has_tm_helix(protein):
      n = len(protein['%s_helices' % program])
      details += [program + "(%d)" % n]

  if is_lipop:
    chop_nterminal_peptide(protein, i_lipop_cut)
  elif is_signalp:
    chop_nterminal_peptide(protein, i_signalp_cut)

  if is_hmm_profile_match:
    category =  "PSE-Cellwall"
  elif has_tm_helix(protein):
    if has_surface_exposed_loop(protein):
      category = "PSE-Membrane"
    else:
      category = "MEMBRANE(non-PSE)"
  else:
    if is_lipop:
      # whole protein considered outer terminal loop
      if sequence_length(protein) < terminal_exposed_loop_min:
        category = "LIPOPROTEIN(non-PSE)"
      else:
        category = "PSE-Lipoprotein"
    elif is_signalp:
      category = "SECRETED"
    else:
      category = "CYTOPLASM(non-PSE)"

  if details == []:
    details = ["."]

  protein['details'] = details
  protein['category'] = category

  return details, category


def protein_output_line(seqid, proteins):
  return '%-15s   %-18s  %-50s  %s' % \
      (seqid, 
      proteins[seqid]['category'], 
      ";".join(proteins[seqid]['details']),
      proteins[seqid]['name'][:60])

def protein_csv_line(seqid, proteins):
  return '%s,%s,%s,"%s"\n' % \
      (seqid, 
       proteins[seqid]['category'], 
       ";".join(proteins[seqid]['details']),
       proteins[seqid]['name'])


def summary_table(params, proteins):
  """
  Returns a string representing a simple summary table of
  protein classifcations.
  """
  out = ""
  counts = {}
  for seqid in proteins:
    category = proteins[seqid]['category']
    
    if category not in counts:
      counts[category] = 1
    else:
      counts[category] += 1
      
  out += "\n\n# Number of proteins in each class:\n"
  pse_total = 0
  
  for c in counts:
    if "PSE-" in c:
      pse_total += counts[c]
  counts["PSE(total)"] = pse_total
  
  for c in sorted(counts):
    if "PSE-" in c:
      spacer = "  "
      pse_total += counts[c]
    else:
      spacer = ""
    out += "# %s%-15s\t%i\n" % (spacer, c, counts[c])
  
  #out += "#\n# PSE(total)\t\t%i\n" % (pse_total)

  return out
