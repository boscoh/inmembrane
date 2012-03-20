from inmembrane import run, parse_fasta_header


def tmhmm(params, proteins):
  tmhmm_out = 'tmhmm.out'
  run('%(tmhmm_bin)s %(fasta)s' % params, tmhmm_out)
  name = None
  for i_line, l in enumerate(open(tmhmm_out)):
    if i_line == 0:
      continue
    words = l.split()
    if not words:
      continue
    if l.startswith("#"):
      name = parse_fasta_header(words[1])[0]
    else:
      name = parse_fasta_header(words[0])[0]
    if name is None:
      continue
    if 'tmhmm_helices' not in proteins[name]:
      proteins[name].update({
        'sequence_length':0,
        'tmhmm_helices':[],
        'tmhmm_inner_loops':[],
        'tmhmm_outer_loops':[]
      })
    if 'Number of predicted TMHs' in l:
      n_helix = int(words[-1])
    if 'Length' in l:
      proteins[name]['sequence_length'] = int(words[-1])
    if 'inside' in l:
      proteins[name]['tmhmm_inner_loops'].append(
          (int(words[-2]), int(words[-1])))
    if 'outside' in l:
      proteins[name]['tmhmm_outer_loops'].append(
          (int(words[-2]), int(words[-1])))
    if 'TMhelix' in l:
      proteins[name]['tmhmm_helices'].append(
          (int(words[-2]), int(words[-1])))
  return proteins

