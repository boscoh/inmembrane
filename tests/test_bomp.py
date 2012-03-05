import os
import unittest
from plugins.bomp_web import *
from inmembrane import create_protein_data_structure, get_params

class TestBomp(unittest.TestCase):
  def setUp(self):
    self.parms = get_params()
    # TODO: put output file in the right place
    self.parms['fasta'] = "tests/bomp/bomps.fasta"
    self.expected_output = {u'gi|107837101': 3, \
                            u'gi|107836588': 5, \
                            u'gi|107836852': 5}
    self.prot_ids, \
    self.proteins = create_protein_data_structure(self.parms['fasta'])
    # run bomp
    self.output = bomp_web(self.parms, self.proteins, force=True)

    print "####   BOMP   #####"
    print self.output
    print "###################"

  def test_bomp(self):
    self.assertEqual(self.expected_output, self.output)

if __name__ == '__main__':
  unittest.main()
