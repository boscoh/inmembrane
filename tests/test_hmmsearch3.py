import os
import unittest
import inmembrane

class TestHmmsearch3(unittest.TestCase):
  def setUp(self):
    self.parms = inmembrane.get_params()
    # TODO: put output file in the right place
    self.parms['fasta'] = "tests/hmmsearch3/hmmsearch3.fasta"
    self.expected_output = {u'SPy_0128': ['LPxTG'], \
                            u'SPy_0191a': ['SLH_ls'], \
                            }
    self.prot_ids, \
    self.proteins = inmembrane.create_protein_data_structure(self.parms['fasta'])
    inmembrane.hmmsearch3(self.parms, self.proteins)

  def test_hmmsearch3(self):
    for seqid in self.expected_output:
        self.assertEqual(self.expected_output[seqid], self.proteins[seqid]['hmmsearch'])

if __name__ == '__main__':
  unittest.main()
