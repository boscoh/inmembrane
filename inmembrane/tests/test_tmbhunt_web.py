import os
import unittest
import sys

import inmembrane
import inmembrane.tests
from inmembrane import helpers
from inmembrane.plugins import tmbhunt_web

from inmembrane.tests.PluginTestBase import PluginTestBase

class TestTmbhunt(PluginTestBase):
  _plugin_name = "tmbhunt_web"
  
  def test_tmbhunt_web(self):
    self.expected_output = {
        'gi|107836852': {
            'tmbhunt_prob': 0.955956, 
            'tmbhunt': True
        }, 
        'gi|107837106': {
            'tmbhunt_prob': 0.23573599999999995, 
            'tmbhunt': False
        }, 
        'gi|107837101': {
            'tmbhunt_prob': 0.955956, 
            'tmbhunt': True
        }, 
        'gi|107836588': {
            'tmbhunt_prob': 0.903904, 
            'tmbhunt': True
        }, 
        'gi|107837107': {
            'tmbhunt_prob': 0.011001000000000039, 
            'tmbhunt': False
        }
    }

    # run TMB-HUNT
    self.output = tmbhunt_web.annotate(self.params, self.proteins, force=True)
    #print self.expected_output
    #helpers.print_proteins(self.output)
    for seqid in self.expected_output:
      self.assertIn(seqid, self.expected_output)
      if seqid in self.output:
        self.assertEqual(self.expected_output[seqid], self.output[seqid])

if __name__ == '__main__':
  unittest.main()
