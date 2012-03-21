from helpers import *


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


def get_annotations(params):
  annotations = [ \
      'annotate_signalp4', 'annotate_lipop1', 
      'annotate_hmmsearch3']

  if dict_get(params, 'helix_programs'):
    if 'tmhmm' in params['helix_programs']:
      annotations.append('annotate_tmhmm')
    if 'memsat3' in params['helix_programs']:
      annotations.append('annotate_memsat3')

  params['hmm_profiles_dir'] = os.path.join(
      os.path.dirname(__file__), 'gram_pos_profiles')

  return annotations


def post_process_protein(params, protein):

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

  if details.endswith(';'):
    details = details[:-1]
  if details is '':
    details = "."

  protein['details'] = details
  protein['category'] = category

  return details, category


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
       proteins[seqid]['name'][:60])







