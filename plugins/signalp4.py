import helpers

def signalp4(params, proteins):
  signalp4_out = 'signalp.out'
  helpers.run('%(signalp4_bin)s -t %(organism)s  %(fasta)s' % params, signalp4_out)
  for l in open(signalp4_out):
    if l.startswith("#"):
      continue
    words = l.split()
    name = helpers.parse_fasta_header(">"+words[0])[0]
    proteins[name].update({ 
      'is_signalp': (words[9] == "Y"),
      'signalp_cleave_position': int(words[4]),
    })
    
