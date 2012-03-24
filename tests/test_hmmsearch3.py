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

    inmembrane.silence_log(True)
    inmembrane.clean_directory('.', ['input.fasta'])
    
    self.params = inmembrane.get_params()
    self.params['fasta'] = "input.fasta"
    self.params['hmm_profiles_dir'] = "../../protocols/gram_pos_profiles"
    self.seqids, self.proteins = \
        inmembrane.create_proteins_dict(self.params['fasta'])
    inmembrane.annotate_hmmsearch3(self.params, self.proteins)

    self.expected_output = {
        u'SPy_0128': ['LPxTG'], 
        u'SPy_0191a': ['SLH'], 
    }
    for seqid in self.expected_output:
      for motif in self.expected_output[seqid]:
        self.assertTrue(motif in self.proteins[seqid]['hmmsearch'])

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
