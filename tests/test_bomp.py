import os
import unittest
from inmembrane import bomp_web, create_protein_data_structure, get_params

class TestBomp(unittest.TestCase):
  def setUp(self):
    self.parms = get_params()
    # TODO: put output file in the right place
    self.parms['fasta'] = "tests/bomp/bomps.fasta"
    self.expected_bomp_categories = {u'107837101': 3, \
                                     u'107836588': 5, \
                                     u'107836852': 5}
    self.prot_ids, \
    self.proteins = create_protein_data_structure(self.parms['fasta'])
    # run bomp
    self.bomp_categories = bomp_web(self.parms, self.proteins, force=True)

  def test_bomp(self):
    self.assertEqual(self.expected_bomp_categories, self.bomp_categories)

if __name__ == '__main__':
  unittest.main()
