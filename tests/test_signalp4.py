import os
import unittest

# hack to allow tests to find inmembrane in directory above
import sys
sys.path.insert(0, '..')

import inmembrane 


class TestSignalp(unittest.TestCase):
  def setUp(self):
    top_tests_dir = os.path.dirname(__file__)
    test_dir = os.path.join(top_tests_dir, 'signalp4')
    self.dir = os.path.abspath(test_dir)

  def test_signalp(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)
    os.system('rm *out')

    self.params = inmembrane.get_params()
    self.params['fasta'] = "signalp4.fasta"
    self.prot_ids, self.proteins = \
        inmembrane.create_protein_data_structure(self.params['fasta'])
    inmembrane.signalp4(self.params, self.proteins)

    self.expected_output = {u'SPy_0252': True, \
                            u'SPy_2077': False, \
                            u'SPy_0317': True}
    for prot_id in self.expected_output:
      self.assertEqual(
          self.expected_output[prot_id], self.proteins[prot_id]['is_signalp'])

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
