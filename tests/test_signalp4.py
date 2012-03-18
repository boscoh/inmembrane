import unittest
import inmembrane 

class TestSignalp(unittest.TestCase):
  def setUp(self):
    self.parms = inmembrane.get_params()
    # TODO: put output file in the right place
    self.parms['fasta'] = "tests/signalp4/signalp4.fasta"
    self.expected_output = {u'SPy_0252': True, \
                            u'SPy_2077': False, \
                            u'SPy_0317': True}
    self.prot_ids, \
    self.proteins = inmembrane.create_protein_data_structure(self.parms['fasta'])

    inmembrane.signalp4(self.parms, self.proteins)


  def test_signalp_signal(self):
    for prot_id in self.expected_output:
      self.assertEqual(
          self.expected_output[prot_id], self.proteins[prot_id]['is_signalp'])


if __name__ == '__main__':
  unittest.main()
