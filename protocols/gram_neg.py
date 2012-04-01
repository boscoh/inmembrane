import os
import helpers
from helpers import dict_get, eval_surface_exposed_loop, \
                    chop_nterminal_peptide, log_stderr

def predict_surface_exposure_barrel(params, protein):
  # TODO: This is a placeholder for a function which will do something
  #       similar to predict_surface_exposure, but focused on inferring 
  #       outer membrane beta barrel topology.
  #
  #       At the moment it's not implemented since I have not
  #       found an easily usable method which gives predictions
  #       accurate enough to be of any real use, particularly
  #       in terms of detecting barrel vs. non-barrel regions
  #       in multidomain BB OMPs (TMBETA-NET strand prediction is 
  #       implemented but it suffers from this problem). Probably 
  #       we will have to do domain detection independently, then 
  #       feed only the barrel domain to the strand predictor. 
  #       I've found PSI-PRED actually works reasonably well on OMP BBs - 
  #       this may be a simpler option over requiring a ProfTMB dependency 
  #       or interfacing with the transFold web service.
  #       
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
      counts[category] = 0
    else:
      counts[category] += 1
      
  out += "\n\n# Number of proteins in each class:\n"
  for c in counts:
    out += "# %-15s %i\n" % (c, counts[c])
    
  return out
    
    
def get_annotations(params):
  params['signalp4_organism'] = 'gram-'
  
  annotations = ['annotate_signalp4', 'annotate_lipop1', \
                 'annotate_tatfind_web']
  annotations += ['annotate_bomp_web']
  
  # TMBETA-NET knows to only run on predicted barrels
  # with the category 'OM(barrel)'
  if 'tmbeta' in dict_get(params, 'barrel_programs'):
    annotations.append('annotate_tmbeta_net_web')

  if dict_get(params, 'helix_programs'):
    if 'tmhmm' in params['helix_programs']:
      annotations.append('annotate_tmhmm')
    if 'memsat3' in params['helix_programs']:
      annotations.append('annotate_memsat3')

  # run some hmm profiles to detect features (eg Tat signal)
  annotations += ['annotate_hmmsearch3']
  params['hmm_profiles_dir'] = os.path.join(
      os.path.dirname(__file__), 'gram_neg_profiles')

  return annotations

def post_process_protein(params, protein, stringent=False):
    
  def has_tm_helix(protein):
    for program in params['helix_programs']:
      if dict_get(protein, '%s_helices' % program):
        return True
    return False

  # we use these functions to detect if and TM-containing IM proteins
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
  
  # in terms of most sublocalization logic, a Tat signal is essentially the 
  # same as a Sec (signalp) signal. We use is_signal_pept to denote that 
  # either is present.
  is_signal_pept = False
  if is_signalp or is_tatfind or \
     ("Tat_PS51318" in dict_get(protein, 'hmmsearch')):
    is_signal_pept = True
  
  is_barrel = False
  if dict_get(protein, 'bomp') >= params['bomp_cutoff']:
    is_barrel = True
    details += ['bomp']
  if dict_get(protein, 'tmbhunt_prob') >= params['tmbhunt_cutoff']:
    is_barrel = True 
    details += ['tmbhunt']
  
  # if stringent, predicted OM barrels must also have a predicted
  # signal sequence
  if stringent and is_signal_pept and is_barrel:
    category = 'OM(barrel)'
  elif is_barrel:
    category = 'OM(barrel)'
  protein['category'] = category

  # set number of predicted OM barrel strands in details
  if (dict_get(protein, 'category') == 'OM(barrel)') and \
      dict_get(protein, 'tmbeta_strands'):
    num_strands = len(protein['tmbeta_strands'])
    details += ['tmbeta_strands(%i)' % (num_strands)]
  
  if is_signal_pept and not is_lipop:
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
    details += ["hmm(%s)" % ",".join(protein['hmmsearch'])]

  if has_tm_helix(protein) and not is_barrel:
    for program in params['helix_programs']:
      n = len(protein['%s_helices' % program])
      details += [program + "(%d);" % n]
    
    category = "IM"
    if long_in_periplasm(protein):
      category += "+peri"
    if long_in_cytoplasm(protein):
      category += "+cyto"
  elif not is_barrel:
    if is_lipop:
      if dict_get(protein, 'lipop_im_retention_signal'):
        category = "LIPOPROTEIN(IM)"
      else:
        category = "LIPOPROTEIN(OM)"
      pass
    elif (is_signal_pept):
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
