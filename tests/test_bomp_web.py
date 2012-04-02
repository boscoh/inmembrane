import os
import glob
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane
import plugins

class TestBomp(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(module_dir, 'bomp')

  def test_bomp(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    inmembrane.silence_log(True)

    inmembrane.clean_directory('.', ['input.fasta'])

    self.params = inmembrane.get_params()
    self.params['fasta'] = "input.fasta"
    self.expected_output = {
        u'gi|107837101': 3, 
        u'gi|107836588': 5, 
        u'gi|107836852': 5
    }
    self.seqids, self.proteins = \
        inmembrane.create_proteins_dict(self.params['fasta'])

    self.output = plugins.bomp_web.annotate(self.params, self.proteins, force=True)


    self.assertEqual(self.expected_output, self.output)

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
