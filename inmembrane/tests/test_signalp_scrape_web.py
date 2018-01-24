import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import signalp_scrape_web

from inmembrane.tests.PluginTestBase import PluginTestBase


class TestSignalpScrapeWeb(PluginTestBase):
    _plugin_name = "signalp_scrape_web"

    def test_signalp4(self):
        if not self.params['signalp4_bin']:
            self.params['signalp4_bin'] = 'signalp_scrape_web'
        self.params['fasta'] = "input.fasta"
        self.params['signalp4_organism'] = 'gram+'
        signalp_scrape_web.annotate(self.params, self.proteins)

        self.expected_output = {
            u'SPy_0252': True,
            u'SPy_2077': False,
            u'SPy_0317': True,
            u'sp|B7LNW7': True,
        }
        for seqid in self.expected_output:
            self.assertEqual(
                self.expected_output[seqid], self.proteins[seqid]['is_signalp'])
        self.assertEqual(self.proteins[u'sp|B7LNW7']['signalp_cleave_position'],
                         22)


if __name__ == '__main__':
    unittest.main()
