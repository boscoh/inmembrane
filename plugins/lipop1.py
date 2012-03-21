import inmembrane 

def annotate_lipop1(params, proteins):
  lipop1_out = 'lipop.out'
  inmembrane.run('%(lipop1_bin)s %(fasta)s' % params, lipop1_out)
  for seqid in proteins:
    proteins[seqid]['is_lipop'] = False
    proteins[seqid]['lipop_cleave_position'] = None
  for l in open(lipop1_out):
    words = l.split()
    if 'SpII score' in l:
      seqid = inmembrane.parse_fasta_header(words[1])[0]
      if 'cleavage' in l:
        pair = words[5].split("=")[1]
        i = int(pair.split('-')[0])
      else:
        i = None
      proteins[seqid]['is_lipop'] = 'Sp' in words[2]
      proteins[seqid]['lipop_cleave_position'] = i
  return proteins
