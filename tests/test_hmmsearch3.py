import os
import unittest
import sys

# hack to allow tests to find inmembrane in directory above
module_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.insert(0, os.path.join(module_dir, '..'))

import inmembrane 
import plugins

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
    plugins.hmmsearch3.annotate(self.params, self.proteins)

    self.expected_output = {
      'SPy_0191a': {
        'sequence_length': 69, 
        'hmmsearch': [], 
        'name': 'SPy_0191a from AE004092', 
        'seq': 'MSKQEVRDSLSTVVGDLSLTTRENQIGSLFLDVQSDEDFGFKVVKVLKSKGIVLNALDESVCGFKFVVE', 
      },
      'SPy_0128': {
        'sequence_length': 340, 
        'hmmsearch': ['LPxTG'], 
        'name': 'SPy_0128 from AE004092', 
        'seq': 'MKLRHLLLTGAALTSFAATTVHGETVVNGAKLTVTKNLDLVNSNALIPNTDFTFKIEPDTTVNEDGNKFKGVALNTPMTKVTYTNSDKGGSNTKTAEFDFSEVTFEKPGVYYYKVTEEKIDKVPGVSYDTTSYTVQVHVLWNEEQQKPVATYIVGYKEGSKVPIQFKNSLDSTTLTVKKKVSGTGGDRSKDFNFGLTLKANQYYKASEKVMIEKTTKGGQAPVQTEASIDQLYHFTLKDGESIKVTNLPVGVDYVVTEDDYKSEKYTTNVEVSPQDGAVKNIAGNSTEQETSTDKDMTITFTNKKDFEVPTGVAMTVAPYIALGIVAVGGALYFVKKKNA', 
      },
    }

    #import helpers
    #helpers.print_proteins(self.proteins)
    
    self.assertEqual(self.proteins, self.expected_output)

    os.chdir(save_dir)


if __name__ == '__main__':
  unittest.main()
