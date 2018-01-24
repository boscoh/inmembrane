import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import tmhmm_scrape_web

import inmembrane
import inmembrane.plugins as plugins

from inmembrane.tests.PluginTestBase import PluginTestBase


class TestTmhmmWeb(PluginTestBase):
    _plugin_name = "tmhmm_scrape_web"

    def test_tmhmm_scrape_web(self):
        tmhmm_scrape_web.annotate(self.params, self.proteins)

        # helpers.print_proteins(self.proteins)

        self.expected_output = {
            'SPy_1392': {
                'tmhmm_outer_loops':
                    [(30, 43), (94, 102), (160, 163), (244, 257), (307, 310),
                     (366, 374)],
                'name': 'SPy_1392 from AE004092',
                'tmhmm_helices':
                    [(7, 29), (44, 66), (73, 93), (103, 125), (137, 159),
                     (164, 186), (221, 243), (258, 280), (287, 306), (311, 330),
                     (343, 365), (375, 394)],
                'tmhmm_inner_loops':
                    [(1, 6), (67, 72), (126, 136), (187, 220), (281, 286),
                     (331, 342), (395, 398)],
            },
            'SPy_1379': {
                'tmhmm_outer_loops':
                    [(56, 58), (106, 178), (231, 244), (307, 320), (373, 381),
                     (434, 437)],
                'name': 'SPy_1379 from AE004092',
                'tmhmm_helices':
                    [(33, 55), (59, 81), (88, 105), (179, 201), (208, 230),
                     (245, 267), (288, 306), (321, 343), (350, 372), (382, 404),
                     (411, 433)],
                'tmhmm_inner_loops':
                    [(1, 32), (82, 87), (202, 207), (268, 287), (344, 349),
                     (405, 410)],
            },
            'SPy_1949': {
                'tmhmm_outer_loops':
                    [(1, 14), (60, 91), (150, 179), (240, 258), (325, 333),
                     (400, 411)],
                'name': 'SPy_1949 from AE004092',
                'tmhmm_helices':
                    [(15, 30), (37, 59), (92, 114), (127, 149), (180, 197),
                     (217, 239), (259, 281), (302, 324), (334, 356),
                     (377, 399)],
                'tmhmm_inner_loops':
                    [(31, 36), (115, 126), (198, 216), (282, 301), (357, 376)],
            }
        }
        for seqid in self.expected_output:
            self.assertTrue(seqid in self.proteins)
            for prop in self.expected_output[seqid]:
                self.assertEqual(
                    self.expected_output[seqid][prop],
                    self.proteins[seqid][prop])


if __name__ == '__main__':
    unittest.main()
