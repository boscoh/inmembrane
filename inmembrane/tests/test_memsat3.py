import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import memsat3

from inmembrane.tests.PluginTestBase import PluginTestBase


# NOTE: since memsat3 is dependent on using a particular BLAST sequence
#       database, if your local version of this database differs to the
#       one used to generate the expected results below, this test will
#       fail. There is no straightforward way to solve this, without
#       packaging a potentially large BLAST database with inmembrane.
class TestMemsat3(PluginTestBase):
    _plugin_name = "memsat3"

    def test_memsat3(self):
        memsat3.annotate(self.params, self.proteins)

        self.expected_output = {
            'SPy_1392': {
                'memsat3_scores': [17.99, 23.11, 27.57, 23.5, 10.29, 9.08,
                                   34.04, 32.19, 28.49, 17.9, 19.18, 25.52],
                'memsat3_helices': [(8, 27), (44, 68), (73, 92), (102, 126),
                                    (132, 150), (153, 172), (220, 243),
                                    (258, 281), (286, 309), (312, 334),
                                    (345, 368), (373, 392)],
                'name': 'SPy_1392 from AE004092',
                'seq': 'MEKTKRYIIATAGILLHLMLGSTYAWSVYRNPILQETGWDQAPVAFAFSLAIFCLGLSAAFMGNLVEQYGPRLTGTVSAILYASGNMLTGLAIDRKEIWLLYIGYGVIGGLGLGAGYITPISTIIKWFPDKRGMATGFAIMGFGFASLLTSPIAQWLIETEGLVATFYLLGLIYLIVMLFASQLIIKPTAAEIAILDKKRLQNNSYLIEGMTAKEALKTKSFYCLWVILFINITCGLGLISVVAPMAQDLTGMSPEMSAIVVGAMGIFNGFGRLVWASLSDYIGRRVTVILLFLVSIIMTISLIFAHSSLIFMISIATLMTCYGAGFSLIPPYLSDLFGAKELATLHGYILTAWAIAALTGPMLLSITVEWTHNYLLTLCVFIVLYILGLMVALRLKK',
                'memsat3_outer_loops': [(28, 43), (93, 101), (151, 152),
                                        (244, 257), (310, 311), (369, 372)],
                'sequence_length': 398,
                'memsat3_inner_loops': [(1, 7), (69, 72), (127, 131),
                                        (173, 219), (282, 285), (335, 344),
                                        (393, 398)],
            },
            'SPy_1379': {
                'memsat3_scores': [15.02, 20.87, 22.47, 28.1],
                'memsat3_helices': [(35, 59), (249, 272), (286, 309),
                                    (323, 346)],
                'name': 'SPy_1379 from AE004092',
                'seq': 'MTIIIMDSNSAHETDNLSVSFLNFCYNSLMKRHFLLLTFYLFLTGLTAGLVAFILTKAIHLIQSLSFGFSQGSFSTMIASVPPQRRALSLLFAGLLAGLGWHLLAKKGKDIQSIQQIIQDDISFSPWTQFWHGWLQLTTVSMGAPVGREGASREVAVTLTSLWSQRCNLSKADQKLLLACASGAALGAVYNAPLATILFILEAILNRWSLKNIYAACLTSYVAVETVALLQGRHEIQYLMPQQHWTLGTLIGSVLAGLILSLFAHAYKHLLKHLPKADAKSQWFIPKVLIAFSLIAGLSIFFPEILGNGKAGLLFFLHEEPHLSYISWLLVAKAVAISLVFASGAKGGKIAPSMMLGGASGLLLAILSQYLIPLSLSNTLAIMVGATIFLGVINKIPLAAPVFLVEITGQSLLMIIPLALANLIFYFSYQFYRFILK',
                'memsat3_outer_loops': [(60, 248), (310, 322)],
                'sequence_length': 437,
                'memsat3_inner_loops': [(1, 34), (273, 285), (347, 437)],
            },
            'SPy_1949': {
                'memsat3_scores': [11.44, 13.03, 19.61, 22.59, 8.81, 12.86,
                                   18.54, 17.77],
                'memsat3_helices': [(16, 35), (40, 59), (89, 112), (135, 159),
                                    (220, 243), (261, 283), (310, 334),
                                    (337, 356)],
                'name': 'SPy_1949 from AE004092',
                'seq': 'MEALLSFIRDILKEPAFLMGLIAFAGLVALKTPAHKVLTGTLGPILGYLMLVAGAGVIVTNLDPLAKLIEHGFSITGVVPNNEAVTSVAQKILGVETMSILVVGLLLNLAFARFTRFKYIFLTGHHSFFMACLLSAVLGAVGFKGSLLIILDGFLLGAWSAISPAIGQQYTLKVTDGDEIAMGHFGSLGYYLSAWVGSKVGKDSKDTEDLQISEKWSFLRNTTISTGLIMVIFYLVATVASVLRNASVAEELAAGQNPFIFAIKSGLTFAVGVAIVYAGVRMILADLIPAFQGIANKLIPNAIPAVDCAVFFPYAPTAVIIGFASSFVGGLLGMLILGVAGGVLIIPGMVPHFFCGATAEIFGNSTGGRRGAMIGASLMAYYSPSCQPCFYLYLVNLVFQTRPLEMWISVF',
                'memsat3_outer_loops': [(1, 15), (60, 88), (160, 219),
                                        (284, 309), (357, 411)],
                'sequence_length': 411,
                'memsat3_inner_loops': [(36, 39), (113, 134), (244, 260),
                                        (335, 336)],
            },
        }

        for seqid in self.expected_output:
            self.assertTrue(seqid in self.proteins)
            for prop in self.expected_output[seqid]:
                self.assertEqual(
                    self.expected_output[seqid][prop],
                    self.proteins[seqid][prop])


if __name__ == '__main__':
    unittest.main()
