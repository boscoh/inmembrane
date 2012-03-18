#!/usr/bin/env python
import os
import glob
import unittest

test_names = [ \
    os.path.basename(f)[:-3] 
    for f in glob.glob('tests/test*.py')]
for test_name in test_names:
  exec('from tests.%s import *' % test_name)
unittest.main()
