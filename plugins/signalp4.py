import helpers

def signalp4(params, proteins):
  for seqid in proteins:
    proteins[seqid]['is_signalp'] = False
    proteins[seqid]['signalp_cleave_position'] = None

  signalp4_out = 'signalp.out'
  helpers.run('%(signalp4_bin)s -t %(organism)s  %(fasta)s' % params, signalp4_out)

  for l in open(signalp4_out):
    if l.startswith("#"):
      continue
    words = l.split()
    seqid = helpers.parse_fasta_header(">"+words[0])[0]
    if (words[9] == "Y"):
      proteins[seqid]['is_signalp'] = True
      proteins[seqid]['signalp_cleave_position'] = int(words[4])

  return proteins
    
