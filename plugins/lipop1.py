import inmembrane 

def lipop1(params, proteins):
  lipop1_out = 'lipop.out'
  inmembrane.run('%(lipop1_bin)s %(fasta)s' % params, lipop1_out)
  for l in open(lipop1_out):
    words = l.split()
    if 'SpII score' in l:
      name = inmembrane.parse_fasta_header(words[1])[0]
      if 'cleavage' in l:
        pair = words[5].split("=")[1]
        i = int(pair.split('-')[0])
      else:
        i = None
      proteins[name].update({
        'is_lipop': 'Sp' in words[2],
        'lipop_cleave_position': i,
      })
  return proteins
