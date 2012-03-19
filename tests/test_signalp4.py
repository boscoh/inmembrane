import os
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane 


class TestSignalp(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(module_dir, 'signalp4')

  def test_signalp(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)
    os.system('rm *out')

    self.params = inmembrane.get_params()
    self.params['fasta'] = "signalp4.fasta"
    self.prot_ids, self.proteins = \
        inmembrane.create_protein_data_structure(self.params['fasta'])
    inmembrane.signalp4(self.params, self.proteins)

    self.expected_output = {
        u'SPy_0252': True, 
        u'SPy_2077': False, 
        u'SPy_0317': True
    }
    for prot_id in self.expected_output:
      self.assertEqual(
          self.expected_output[prot_id], self.proteins[prot_id]['is_signalp'])

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
