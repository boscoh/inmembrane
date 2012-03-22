import os
import helpers
from helpers import dict_get, eval_surface_exposed_loop, \
                    chop_nterminal_peptide, log_stderr

# the formatting of the annotations printed to stdout and written
# to csv file can be changed by overriding these functions if desired
protein_output_line = helpers.protein_output_line
protein_csv_line = helpers.protein_csv_line

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
