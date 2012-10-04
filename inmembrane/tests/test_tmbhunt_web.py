import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import tmbhunt_web

class TestTmbhunt(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(
       os.path.abspath(
       os.path.dirname(inmembrane.tests.__file__)), 'tmbhunt')

  def test_tmbhunt(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    helpers.silence_log(True)
    helpers.clean_directory('.', ['input.fasta'])
     
    self.params = inmembrane.get_params()
    self.params['fasta'] = "input.fasta"
    self.expected_output = {
        'gi|107836852': {
            'tmbhunt_prob': 0.955956, 
            'tmbhunt': True
        }, 
        'gi|107837106': {
            'tmbhunt_prob': 0.23573599999999995, 
            'tmbhunt': False
        }, 
        'gi|107837101': {
            'tmbhunt_prob': 0.955956, 
            'tmbhunt': True
        }, 
        'gi|107836588': {
            'tmbhunt_prob': 0.903904, 
            'tmbhunt': True
        }, 
        'gi|107837107': {
            'tmbhunt_prob': 0.011001000000000039, 
            'tmbhunt': False
        }
    }
    self.seqids, \
    self.proteins = helpers.create_proteins_dict(self.params['fasta'])

    # run TMB-HUNT
    self.output = tmbhunt_web.annotate(self.params, self.proteins, force=True)
    #print self.expected_output
    #helpers.print_proteins(self.output)
    for seqid in self.expected_output:
      self.assertIn(seqid, self.expected_output)
      if seqid in self.output:
        self.assertEqual(self.expected_output[seqid], self.output[seqid])

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
