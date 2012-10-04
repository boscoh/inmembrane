#!/usr/bin/env python
# Copyright (c) 2012, Bosco Ho <boscoh@gmail.com>
# Copyright (c) 2012, Andrew Perry <ajperry@pansapiens.com>
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or
# without modification, are permitted provided that the following 
# conditions are met:

# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the 
# documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS 
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE 
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, 
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS 
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED 
# AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF 
# THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
# DAMAGE.

import sys
import os
import glob

import inmembrane
#from inmembrane import helpers
from inmembrane.helpers import *
# will load all plugins in the plugins/ directory
from inmembrane.plugins import *

description = """
inmembrane %s (https://github.com/boscoh/inmembrane)

Inmembrane is a proteome annotation pipeline. It takes 
a FASTA file, then carries out sequential analysis of 
each sequence with a bunch of third-party programs, and 
collates the results.

(c) 2011 Bosco Ho and Andrew Perry
""" % (inmembrane.__version__)

if __name__ == "__main__":
  from optparse import OptionParser
  usage = "usage: %prog [options] [sequences.fasta]"
  parser = OptionParser(usage=usage)
  parser.add_option("-t", "--test",
                     action="store_true", dest="run_tests", default=False,
                    help="Run tests checking the each plugin works as expected.\
                          No sequence file should be specified if this flag is\
                          used.")
  parser.add_option("-n", "--tests-no-network",
                     action="store_true", dest="no_network", default=False,
                     help="don't run tests that require a network connection")                   
  (options, args) = parser.parse_args()
  
  if not options.run_tests:
    ##
    ## Setup params and run the analysis
    ## The high level plugin loop is in inmembrane/__init__.py
    ##  
    log_stderr("inmembrane %s (https://github.com/boscoh/inmembrane)"
              % (inmembrane.__version__))
    params = inmembrane.get_params()
    if ('fasta' not in params or not params['fasta']) and not args:
      print description
      parser.print_help()
      sys.exit(1)
    if 'fasta' not in params or not params['fasta']:
      params['fasta'] = args[0]
    inmembrane.process(params)
    
  else:
    ##
    ## Run unit tests instead of 'normal' analysis
    ##
    if len(args) > 0:
      print description
      parser.print_help()
      sys.exit(1)
    # unusual hack to remove our custom flags, since 
    # unittest expects it's own set
    if "-t" in sys.argv:
      sys.argv.remove("-t")
    if "--test" in sys.argv:
      sys.argv.remove("--test")   
    if "-n" in sys.argv:
      sys.argv.remove("-n")
    if "--tests-no-network" in sys.argv:
      sys.argv.remove("--tests-no-network")
      
    file_tag = os.path.join(inmembrane.module_dir, 'tests', 'test*.py')
    test_names = [ \
        os.path.basename(f)[:-3] for f in glob.glob(file_tag)]

    tests_to_run = []
    for test_name in test_names:
      if not (options.no_network and (test_name[-4:] == "_web")):
        exec('from inmembrane.tests.%s import *' % test_name)
        tests_to_run.append(test_name)

    log_stderr("Will run: " + ", ".join(tests_to_run))

    unittest.main()

