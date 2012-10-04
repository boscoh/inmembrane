# -*- coding: utf-8 -*-

citation = {'ref': u"Gromiha MM, Ahmad S, Suwa M. Neural network-based "
                   u"prediction of transmembrane beta-strand segments in "
                   u"outer membrane proteins. Journal of computational "
                   u"chemistry 2004 Apr;25(5):762-7. \n"
                   u"<http://dx.doi.org/10.1002/jcc.10386>",
            'name': "TMBETA-NET"
           }

__DEBUG__ = False

import sys, os, time, StringIO
import json
import twill
from twill.commands import find, formfile, follow, fv, go, show, \
                             showforms, showlinks, submit, agent
                             
from inmembrane.helpers import dict_get, log_stderr

def annotate(params, proteins, \
                   url="http://psfs.cbrc.jp/tmbeta-net/", \
                   category='OM(barrel)',
                   force=False):
  """
  Uses the TMBETA-NET web service (http://psfs.cbrc.jp/tmbeta-net/) to
  predict strands of outer membrane beta-barrels.
  
  By default, category='BARREL' means prediction will only be run
  on proteins in the set with this category property. To process all
  proteins, change category to None.

  These keys are added to the proteins dictionary: 
    'tmbeta_strands' - a list of lists with paired start and end 
                       residues of each predicted strand. 
                       (eg [[3,9],[14,21], ..etc ])
  """

  # set the user-agent so web services can block us if they want ... :/
  python_version = sys.version.split()[0]
  agent("Python-urllib/%s (twill; inmembrane)" % python_version)
  
  outfile = 'tmbeta_net.out'
  log_stderr("# TMBETA-NET(web) %s > %s" % (params['fasta'], outfile))
  
  tmbeta_strands = {}
  if not force and os.path.isfile(outfile):
    log_stderr("# -> skipped: %s already exists" % outfile)
    fh = open(outfile, 'r')
    tmbeta_strands = json.loads(fh.read())
    fh.close()    
    for seqid in tmbeta_strands:
      proteins[seqid]['tmbeta_strands'] = tmbeta_strands[seqid]
      
    return tmbeta_strands

  # dump extraneous output into this blackhole so we don't see it
  if not __DEBUG__: twill.set_output(StringIO.StringIO())

  for seqid in proteins:
    
    # only run on sequences which match the category filter
    if force or \
       (category == None) or \
       (dict_get(proteins[seqid], 'category') == category):
      pass
    else:
      continue
      
    go(url)
    if __DEBUG__: showforms()
    fv("1","sequence",proteins[seqid]['seq'])
    submit()
    log_stderr("# TMBETA-NET: Predicting strands for %s - %s\n" \
                      % (seqid, proteins[seqid]['name']))
    out = show()
    time.sleep(1)
    
    # parse the web page returned, extract strand boundaries
    proteins[seqid]['tmbeta_strands'] = []
    for l in out.split('\n'):
      if __DEBUG__: log_stderr("## " + l)

      if "<BR>Segment " in l:
        i,j = l.split(":")[1].split("to")
        i = int(i.strip()[1:])
        j = int(j.strip()[1:])
        proteins[seqid]['tmbeta_strands'].append([i,j])

        if __DEBUG__: log_stderr("# TMBETA-NET(web) segments: %s, %s" % (i, j))

    tmbeta_strands[seqid] = proteins[seqid]['tmbeta_strands']

  # we store the parsed strand boundaries in JSON format
  fh = open(outfile, 'w')
  fh.write(json.dumps(tmbeta_strands, separators=(',',':\n')))
  fh.close()

  return tmbeta_strands
