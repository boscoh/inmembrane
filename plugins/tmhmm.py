

import helpers


def annotate_tmhmm(params, proteins):
  """
  Runs THMHMM and parses the output files. Takes a standard 'inmembrane'
  params dictionary and a global proteins dictionary which it populates with
  results.
  
  In the current implementation, this function extracts and feeds sequences to tmhmm
  one by one via a temporary file.
  
  These keys are added to the proteins dictionary: 
    - 'tmhmm_helices', a list of tuples describing the first and last residue
      number of each transmembrane segment; 
    
    - 'tmhmm_scores', a list of confidence scores (floats) for each predicted 
      tm segment;
    
    - 'tmhmm_inner_loops', a list of tuples describing the first and last residue
       number of each predicted internal loop segment;
    
    - 'tmhmm_outer_loops', a list of tuples describing the first and last residue
       number of each predicted outer loop segment;
  """

  tmhmm_out = 'tmhmm.out'
  helpers.run('%(tmhmm_bin)s %(fasta)s' % params, tmhmm_out)

  seqid = None
  for i_line, l in enumerate(open(tmhmm_out)):
    if i_line == 0:
      continue
    words = l.split()
    if not words:
      continue

    if l.startswith("#"):
      seqid = helpers.parse_fasta_header(words[1])[0]
    else:
      seqid = helpers.parse_fasta_header(words[0])[0]
    if seqid is None:
      continue

    # initialize fields in proteins[seqid]
    if 'tmhmm_helices' not in proteins[seqid]:
      proteins[seqid].update({
        'sequence_length':0,
        'tmhmm_helices':[],
        'tmhmm_inner_loops':[],
        'tmhmm_outer_loops':[]
      })

    if 'Length' in l:
      proteins[seqid]['sequence_length'] = int(words[-1])
    if 'inside' in l:
      proteins[seqid]['tmhmm_inner_loops'].append(
          (int(words[-2]), int(words[-1])))
    if 'outside' in l:
      proteins[seqid]['tmhmm_outer_loops'].append(
          (int(words[-2]), int(words[-1])))
    if 'TMhelix' in l:
      proteins[seqid]['tmhmm_helices'].append(
          (int(words[-2]), int(words[-1])))

  return proteins

