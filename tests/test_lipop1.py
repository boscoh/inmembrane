import os
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane 


class TestLipoP(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(module_dir, 'lipop1')

  def test_lipop(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    inmembrane.silence_log(True)
    inmembrane.clean_directory('.', ['input.fasta'])
   
    self.params = inmembrane.get_params()
    self.params['fasta'] = "input.fasta"
    self.seqids, self.proteins = \
        inmembrane.create_protein_data_structure(self.params['fasta'])

    inmembrane.lipop1(self.params, self.proteins)

    self.expected_output = {
        u'SPy_0252': True,
        u'SPy_2077': False, 
        u'SPy_0317': True
    }
    for seqid in self.expected_output:
      self.assertEqual(
          self.expected_output[seqid], self.proteins[seqid]['is_lipop'])

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
