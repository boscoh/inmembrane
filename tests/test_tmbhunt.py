import os
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane 


class TestTmbhunt(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(module_dir, 'tmbhunt')

  def test_tmbhunt(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    self.params = inmembrane.get_params()
    self.params['fasta'] = "bomps.fasta"
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
    self.prot_ids, \
    self.proteins = inmembrane.create_protein_data_structure(self.params['fasta'])

    # run TMB-HUNT
    self.output = inmembrane.tmbhunt_web(self.params, self.proteins, force=True)
    
    self.assertEqual(self.expected_output, self.output)

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
