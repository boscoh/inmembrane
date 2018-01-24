import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import lipop_scrape_web

from inmembrane.tests.PluginTestBase import PluginTestBase


class TestLipoP(PluginTestBase):
    _plugin_name = "lipop_scrape_web"

    def test_lipop_scrape_web(self):
        lipop_scrape_web.annotate(self.params, self.proteins)

        self.expected_output = {
            u'SPy_0252': True,
            u'SPy_2077': False,
            u'SPy_0317': True,
            u'tr|Q9HYX8': True,
        }

        # print self.expected_output
        # print
        # print self.proteins

        for seqid in self.expected_output:
            self.assertEqual(
                self.expected_output[seqid], self.proteins[seqid]['is_lipop'])
        self.assertEqual(self.proteins[u'tr|Q9HYX8']['lipop_cleave_position'],
                         19)
        self.assertIn('lipop_im_retention_signal', self.proteins[u'tr|Q9HYX8'])
        self.assertTrue(
            self.proteins[u'tr|Q9HYX8']['lipop_im_retention_signal'])


if __name__ == '__main__':
    unittest.main()
