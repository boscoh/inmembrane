import os
import unittest
from inmembrane import tmbhunt_web, create_protein_data_structure, get_params

class TestTmbhunt(unittest.TestCase):
  def setUp(self):
    self.parms = get_params()
    # TODO: put output file in the right place
    self.parms['fasta'] = "tests/tmbhunt/bomps.fasta"
    self.expected_tmbhunt_output = {'gi|107837107': {'tmbhunt_prob': 0.988999, 'tmbhunt': False}, 'gi|107837101': {'tmbhunt_prob': 0.955956, 'tmbhunt': True}, 'gi|107836588': {'tmbhunt_prob': 0.903904, 'tmbhunt': True}, 'gi|107837106': {'tmbhunt_prob': 0.764264, 'tmbhunt': False}, 'gi|107836852': {'tmbhunt_prob': 0.955956, 'tmbhunt': True}}
    self.prot_ids, \
    self.proteins = create_protein_data_structure(self.parms['fasta'])
    # run TMB-HUNT
    self.tmbhunt_output = tmbhunt_web(self.parms, self.proteins, force=True)

  def test_tmbhunt(self):
    self.assertEqual(self.expected_tmbhunt_output, self.tmbhunt_output)

if __name__ == '__main__':
  unittest.main()
