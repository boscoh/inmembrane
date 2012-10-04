import os
from inmembrane.helpers import dict_get, eval_surface_exposed_loop, \
                    chop_nterminal_peptide, log_stderr

    
def get_annotations(params):
  annotations = []
  params['signalp4_organism'] = 'gram-'
  
  if not params['signalp4_bin'] or params['signalp4_bin'] == 'signalp_web':
    annotations += ['signalp_web']
  else:
    annotations += ['signalp4']
  
  if not params['lipop1_bin'] or params['lipop1_bin'] == 'lipop_web':
    annotations += ['lipop_web']
  else:
    annotations += ['lipop1']
  
  annotations += ['tatfind_web']
  
  if 'bomp' in dict_get(params, 'barrel_programs'):
    annotations.append('bomp_web')
  if 'tmbhunt' in dict_get(params, 'barrel_programs'):
    annotations.append('tmbhunt_web')
  if 'tmbetadisc-rbf' in dict_get(params, 'barrel_programs'):
    annotations.append('tmbetadisc_rbf_web')
    
  # TMBETA-NET knows to only run on predicted barrels
  # with the category 'OM(barrel)'
  if 'tmbeta' in dict_get(params, 'barrel_programs'):
    annotations.append('tmbeta_net_web')

  if dict_get(params, 'helix_programs'):
    if 'tmhmm' in params['helix_programs']:
      annotations.append('tmhmm')
    if 'memsat3' in params['helix_programs']:
      annotations.append('memsat3')

  # run some hmm profiles to detect features (eg Tat signal)
  annotations += ['hmmsearch3']
  params['hmm_profiles_dir'] = os.path.join(
      os.path.dirname(__file__), 'gram_neg_profiles')

  return annotations

def post_process_protein(params, protein):
    
  def has_tm_helix(protein):
    for program in params['helix_programs']:
      if dict_get(protein, '%s_helices' % program):
        return True
    return False

  # these functions detect if and TM-containing IM proteins
  # have large loops / terminal regions in the periplasm or cytoplasm
  # that may be accessible / inaccessible in spheroplast shaving 
  # experiments.
  def has_long_loops(protein, loop_str='_outer_loops', \
                     loop_length=params['internal_exposed_loop_min']):
    for annot in protein:
      if loop_str in annot:
        for loop in protein[annot]:
          l_len = loop[1]-loop[0]
          if l_len >= loop_length:
            return True
    return False
  
  def long_in_periplasm(protein, \
                        loop_length=params['internal_exposed_loop_min']):
    return has_long_loops(protein, '_outer_loops', loop_length)
  
  def long_in_cytoplasm(protein, \
                        loop_length=params['internal_exposed_loop_min']):
    return has_long_loops(protein, '_inner_loops', loop_length)

  
  details = []
  category = "UNKNOWN"
  is_hmm_profile_match = dict_get(protein, 'hmmsearch')
  is_signalp = dict_get(protein, 'is_signalp')
  is_tatfind = dict_get(protein, 'is_tatfind')
  is_lipop = dict_get(protein, 'is_lipop')
  
  # in terms of most sublocalization logic, a Tat signal is similar to a 
  # Sec (signalp) signal. We use has_signal_pept to denote that either 
  # is present.
  has_signal_pept = False
  if is_signalp or is_tatfind or \
     (('hmmsearch' in protein) and "Tat_PS51318" in protein['hmmsearch']):
    has_signal_pept = True
  
  # annotate the barrels - high scoring bomp hits don't require a 
  # signal peptide, low scoring ones do
  has_barrel = False
  bomp_score = dict_get(protein, 'bomp')
  if (bomp_score >= params['bomp_clearly_cutoff']) or \
     (has_signal_pept and bomp_score >= params['bomp_maybe_cutoff']):
    
    details += ['bomp(%i)' % (bomp_score)]
    has_barrel = True
    
  tmbhunt_prob = dict_get(protein, 'tmbhunt_prob')
  if (tmbhunt_prob >= params['tmbhunt_clearly_cutoff']) or \
     (has_signal_pept and tmbhunt_prob >= params['tmbhunt_maybe_cutoff']):
    details += ['tmbhunt(%.2f)' % (tmbhunt_prob)]
    has_barrel = True
    
  if has_signal_pept and dict_get(protein, 'is_tmbetadisc_rbf'):
    details += ['tmbetadisc-rbf']
    has_barrel = True
    
  if has_barrel:
    category = 'OM(barrel)'
    
  # we only regard the barrel prediction as a true positive
  # if a signal peptide is also present
#  is_barrel = False
#  if has_signal_pept and has_barrel: # TODO and num_tms <= 1:
#    category = 'OM(barrel)'
#    is_barrel = True
    
  # set number of predicted OM barrel strands in details
  if has_barrel and \
      dict_get(protein, 'tmbeta_strands'):
    num_strands = len(protein['tmbeta_strands'])
    details += ['tmbeta_strands(%i)' % (num_strands)]
  
  if has_signal_pept and not is_lipop and \
    (dict_get(protein, 'signalp_cleave_position')):
    # we use the SignalP signal peptidase cleavage site for Tat signals
    chop_nterminal_peptide(protein,  protein['signalp_cleave_position'])
  
  if is_tatfind:
    details += ["tatfind"]
  
  if is_signalp:
    details += ["signalp"]
  
  if is_lipop:
    details += ["lipop"]
    chop_nterminal_peptide(protein, protein['lipop_cleave_position'])
  
  if is_hmm_profile_match:
    details += ["hmm(%s)" % "|".join(protein['hmmsearch'])]

  if has_tm_helix(protein) and not has_barrel:
    for program in params['helix_programs']:
      n = len(protein['%s_helices' % program])
      details += [program + "(%d)" % n]
    
    category = "IM"
    if long_in_periplasm(protein):
      category += "+peri"
    if long_in_cytoplasm(protein):
      category += "+cyto"
  elif not has_barrel:
    if is_lipop:
      if dict_get(protein, 'lipop_im_retention_signal'):
        category = "LIPOPROTEIN(IM)"
      else:
        category = "LIPOPROTEIN(OM)"
      pass
    elif (has_signal_pept):
      category = "PERIPLASMIC/SECRETED"
    else:
      category = "CYTOPLASM"

  if details == []:
    details = ["."]

  protein['details'] = details
  protein['category'] = category

  return details, category

def protein_output_line(seqid, proteins):
  return '%-15s   %-13s  %-50s  %s' % \
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
  for c in counts:
    out += "# %-15s %i\n" % (c, counts[c])
    
  return out
