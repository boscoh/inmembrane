import os
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane 


class TestHmmsearch3(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(module_dir, 'hmmsearch3')

  def test_hmmsearch3(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    self.params = inmembrane.get_params()
    self.params['fasta'] = "hmmsearch3.fasta"
    self.seqids, self.proteins = \
        inmembrane.create_protein_data_structure(self.params['fasta'])
    inmembrane.hmmsearch3(self.params, self.proteins)

    self.expected_output = {
        u'SPy_0128': ['LPxTG'], 
        u'SPy_0191a': ['SLH_ls'], 
    }
    for seqid in self.expected_output:
      for motif in self.expected_output[seqid]:
        self.assertTrue(motif in self.proteins[seqid]['hmmsearch'])

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
