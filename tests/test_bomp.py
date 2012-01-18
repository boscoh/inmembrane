import unittest
from inmembrane import bomp_web, init, parms

class TestBomp(unittest.TestCase):
  def setUp(self):
    self.parms = parms
    self.parms['fasta'] = "tests/bomp/bomps.fasta"
    self.expected_bomp_categories = {u'gi|107837101|gb|ABF84970.1|': 3, u'gi|107836588|gb|ABF84457.1|': 5, u'gi|107836852|gb|ABF84721.1|': 5}
    self.prot_ids, self.proteins = init(parms)
    # run bomp
    self.bomp_categories = bomp_web(self.parms, self.proteins, force=True)

  def test_bomp(self):
    self.assertEqual(self.expected_bomp_categories, self.bomp_categories)

if __name__ == '__main__':
  unittest.main()
