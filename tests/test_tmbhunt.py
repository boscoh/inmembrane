import os
import unittest
from plugins.tmbhunt_web import *
from inmembrane import create_protein_data_structure, get_params

class TestTmbhunt(unittest.TestCase):
  def setUp(self):
    self.parms = get_params()
    # TODO: put output file in the right place
    self.parms['fasta'] = "tests/tmbhunt/bomps.fasta"
    self.expected_output = {'gi|107836852': {'tmbhunt_prob': 0.955956, 'tmbhunt': True}, 'gi|107837106': {'tmbhunt_prob': 0.23573599999999995, 'tmbhunt': False}, 'gi|107837101': {'tmbhunt_prob': 0.955956, 'tmbhunt': True}, 'gi|107836588': {'tmbhunt_prob': 0.903904, 'tmbhunt': True}, 'gi|107837107': {'tmbhunt_prob': 0.011001000000000039, 'tmbhunt': False}}
    self.prot_ids, \
    self.proteins = create_protein_data_structure(self.parms['fasta'])
    # run TMB-HUNT
    self.output = tmbhunt_web(self.parms, self.proteins, force=True)
    
    print "#### TMB-HUNT #####"
    print self.output
    print "###################"
    
  def test_tmbhunt(self):
    self.assertEqual(self.expected_output, self.output)

if __name__ == '__main__':
  unittest.main()
