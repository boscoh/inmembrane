from inmembrane import run, parse_fasta_header

def signalp4(params, proteins):
  signalp4_out = 'signalp.out'
  run('%(signalp4_bin)s -t %(organism)s  %(fasta)s' % params, signalp4_out)
  for l in open(signalp4_out):
    if l.startswith("#"):
      continue
    words = l.split()
    name = parse_fasta_header(">"+words[0])[0]
    proteins[name].update({ 
      'is_signalp': (words[9] == "Y"),
      'signalp_cleave_position': int(words[4]),
    })
    
if __name__ == "__main__":
  import unittest
  from inmembrane.tests import test_signalp
  unittest.main()
