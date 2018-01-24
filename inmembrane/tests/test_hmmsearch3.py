import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import hmmsearch3

from inmembrane.tests.PluginTestBase import PluginTestBase


class TestHmmsearch3(PluginTestBase):
    _plugin_name = "hmmsearch3"

    def test_hmmsearch3(self):
        self.params['hmm_profiles_dir'] = os.path.join(inmembrane.module_dir,
                                                       "protocols/gram_pos_profiles")
        hmmsearch3.annotate(self.params, self.proteins)

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

        # helpers.print_proteins(self.proteins)
        for seqid in self.proteins:
            self.assertEqual(self.proteins[seqid]['hmmsearch'],
                             self.expected_output[seqid]['hmmsearch'])


if __name__ == '__main__':
    unittest.main()
