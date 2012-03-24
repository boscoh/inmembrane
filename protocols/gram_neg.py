import os
import helpers
from helpers import dict_get, eval_surface_exposed_loop, \
                    chop_nterminal_peptide, log_stderr

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


def print_summary_table(params, proteins):
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
      
  log_stderr("# Number of proteins in each class:")
  for c in counts:
    log_stderr("%-15s %i" % (c, counts[c]))

# TODO: Workflow:
#       These rules can be approached in two ways. Run every annotation plugin
#       on every sequence, then post-hoc apply these rules (potentially inefficient
#       and rude to web services). Alternatively the protocol can control the
#       annotation loop and only run certain analyses (eg BOMP) if a certain 
#       condition is met (eg is_signalp). 
#       * If is_lipop, 
#         check retention signal annotation (+2Asp) -> periplasmic inner or outer
#       * If is_signalp (or TAT) and no other TM and not is_lipop -> periplasmic
#       * If is_signalp, run OMPBB (BOMP) check. If OMPBB -> outer membrane bb
#       * If is_signalp, and has after the signal tranmembranes -> inner membrane
#       * If inner membrane, check for long cyto or peri loops: inner[+cyto][+peri]
#       * If not is_signalp and not is_lipop (and not TAT): -> cytoplasmic
def get_annotations(params):
  annotations = ['annotate_signalp4', 'annotate_lipop1']
  #annotations += ['annotate_tatp']
  annotations += ['annotate_bomp']

  if dict_get(params, 'helix_programs'):
    if 'tmhmm' in params['helix_programs']:
      annotations.append('annotate_tmhmm')
    if 'memsat3' in params['helix_programs']:
      annotations.append('annotate_memsat3')

  # currently we don't use any HMM profiles on gram-
  # annotation, but this may be useful in the future
  #annotations += ['annotate_hmmsearch3']
  #params['hmm_profiles_dir'] = os.path.join(
  #    os.path.dirname(__file__), 'gram_neg_profiles')

  return annotations

def post_process_protein(params, protein):
    
  def has_tm_helix(protein):
    for program in params['helix_programs']:
      if dict_get(protein, '%s_helices' % program):
        return True
    return False

  # TODO: rather than this, lets classify each inner membrane protein as 
  # IM with optional large periplasmic or large cytoplasmic loops/domain
  # (eg possible classifications: IM+cyt, IM+peri, IM+cyt+peri )
  def has_periplasmic_loop(protein):
    for program in params['helix_programs']:
      if eval_surface_exposed_loop(
          protein['sequence_length'], 
          len(protein['%s_helices' % (program)]), 
          protein['%s_outer_loops' % (program)], 
          params['terminal_exposed_loop_min'], 
          params['internal_exposed_loop_min']):
        return True
    return False
      
  #is_hmm_profile_match = dict_get(protein, 'hmmsearch')
  is_lipop = dict_get(protein, 'is_lipop')
  if is_lipop:
    i_lipop_cut = protein['lipop_cleave_position']
  is_signalp = dict_get(protein, 'is_signalp')
  if is_signalp:
    i_signalp_cut = protein['signalp_cleave_position']

  # TODO: make details a plain list.
  #       leave the insertion of ';' to output time
  details = ""
  #if is_hmm_profile_match:
  #  details += "hmm(%s);" % protein['hmmsearch'][0]
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

  elif has_tm_helix(protein):
    if has_periplasmic_loop(protein):
      category = "PSE"
    else:
      category = "MEMBRANE"
  else:
    if is_lipop:
      # whole protein considered outer terminal loop
      if len(protein['seq']) < params['terminal_exposed_loop_min']:
        category = "MEMBRANE"
      else:
        category = "PSE"
    elif is_signalp:
      category = "SECRETED"
    else:
      category = "CYTOPLASM"

  if details.endswith(';'):
    details = details[:-1]
  if details is '':
    details = "."

  protein['details'] = details
  protein['category'] = category

  return details, category

#def identify_omps(params, stringent=False):
#  """
#  Identifies outer membrane proteins from gram-negative bacteria.
#  
#  If stringent=True, all predicted outer membrane barrels must also
#  have a predicted signal sequence to be categorized as BARREL.
#  """
#  
#  seqids, proteins = create_proteins_dict(params['fasta'])
#
#  features = [signalp4, lipop1, hmmsearch3]
#  if dict_get(params, 'helix_programs'):
#    if 'tmhmm' in params['helix_programs']:
#      features.append(tmhmm)
#    if 'memsat3' in params['helix_programs']:
#      features.append(memsat3)
#  if dict_get(params, 'barrel_programs'):
#    if 'tmbhunt' in params['barrel_programs']:
#      features.append(annotate_tmbhunt_web)
#    if 'bomp' in params['barrel_programs']:
#      features.append(annotate_bomp_web)
#  for extract_protein_feature in features:
#    extract_protein_feature(params, proteins)
#  
#  for seqid, protein in proteins.items():
#    # TODO: this is used for setting 'category', however
#    #       we may need to make a gram- OM specific version
#    #       (eg, run after strand prediction so we can look at
#    #            strand topology, detect long extracellular loops etc) 
#    details, category = predict_surface_exposure(params, protein)
#    proteins[seqid]['category'] = category
#    proteins[seqid]['details'] = details
#    
#    if stringent:
#      if dict_get(protein, 'is_signalp') and \
#       ( dict_get(protein, 'bomp') or \
#         dict_get(protein, 'tmbhunt') ):
#       proteins[seqid]['category'] = 'BARREL'
#    else:
#      if dict_get(protein, 'bomp') or \
#         dict_get(protein, 'tmbhunt'):
#         proteins[seqid]['category'] = 'BARREL'
#    
#  # TMBETA-NET knows to only run on predicted barrels
#  if 'tmbeta' in params['barrel_programs']:
#    annotate_tmbeta_net_web((params, proteins, category='BARREL')
#
#  for seqid in proteins:
#    details = proteins[seqid]['details']
#    if dict_get(proteins[seqid], 'tmbeta_strands'):
#      num_strands = len(proteins[seqid]['tmbeta_strands'])
#      details += 'tmbeta(%i)' % (num_strands)
#    if details.endswith(';'):
#      details = details[:-1]
#    if details is '':
#      details = "."
#    proteins[seqid]['details'] = details
#    
#  print_summary_table(proteins)
#
#  return seqids, proteins


def protein_output_line(seqid, proteins):
  return '%-15s   %-13s  %-50s  %s' % \
      (seqid, 
      proteins[seqid]['category'], 
      proteins[seqid]['details'],
      proteins[seqid]['name'][:60])

def protein_csv_line(seqid, proteins):
  return '%s,%s,%s,"%s"\n' % \
      (seqid, 
       proteins[seqid]['category'], 
       proteins[seqid]['details'],
       proteins[seqid]['name'])
