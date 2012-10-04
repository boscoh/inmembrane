import os
import glob
import unittest
import sys

import inmembrane
import inmembrane.tests

from inmembrane import helpers
from inmembrane.plugins import bomp_web

class TestBomp(unittest.TestCase):
  def setUp(self):
    self.dir = os.path.join(
               os.path.abspath(
               os.path.dirname(inmembrane.tests.__file__)), 'bomp')

  def test_bomp(self):
    save_dir = os.getcwd()
    os.chdir(self.dir)

    helpers.silence_log(True)

    helpers.clean_directory('.', ['input.fasta'])

    self.params = inmembrane.get_params()
    self.params['fasta'] = "input.fasta"
    self.expected_output = {
        u'gi|107837101': 3, 
        u'gi|107836588': 5, 
        u'gi|107836852': 5
    }
    self.seqids, self.proteins = \
        helpers.create_proteins_dict(self.params['fasta'])

    self.output = bomp_web.annotate(self.params, self.proteins, force=True)


    self.assertEqual(self.expected_output, self.output)

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
