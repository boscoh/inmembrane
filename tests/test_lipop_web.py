import os
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane 
import plugins

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
        inmembrane.create_proteins_dict(self.params['fasta'])

    plugins.lipop_web.annotate(self.params, self.proteins)

    self.expected_output = {
        u'SPy_0252': True,
        u'SPy_2077': False, 
        u'SPy_0317': True,
        u'tr|Q9HYX8' : True,
    }
    
    for seqid in self.expected_output:
      self.assertEqual(
          self.expected_output[seqid], self.proteins[seqid]['is_lipop'])
    self.assertIn('lipop_im_retention_signal', self.proteins[u'tr|Q9HYX8'])
    self.assertTrue(self.proteins[u'tr|Q9HYX8']['lipop_im_retention_signal'])
    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
