import sys, os, time, StringIO
import json
import twill
from twill.commands import find, formfile, follow, fv, go, show, \
                             showforms, showlinks, submit, agent
                             
import inmembrane



def annotate_tmbeta_net_web(params, proteins, \
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
  inmembrane.log_stderr("# TMBETA-NET(web) %s > %s" % (params['fasta'], outfile))
  
  tmbeta_strands = {}
  if not force and os.path.isfile(outfile):
    inmembrane.log_stderr("# -> skipped: %s already exists" % outfile)
    fh = open(outfile, 'r')
    tmbeta_strands = json.loads(fh.read())
    fh.close()    
    for seqid in tmbeta_strands:
      proteins[seqid]['tmbeta_strands'] = tmbeta_strands[seqid]
      
    return tmbeta_strands

  # dump extraneous output into this blackhole so we don't see it
  if not inmembrane.LOG_DEBUG: twill.set_output(StringIO.StringIO())

  for seqid in proteins:
    
    # only run on sequences which match the category filter
    if force or \
       (category == None) or \
       (inmembrane.dict_get(proteins[seqid], 'category') == category):
      pass
    else:
      continue
      
    go(url)
    if inmembrane.LOG_DEBUG: showforms()
    fv("1","sequence",proteins[seqid]['seq'])
    submit()
    inmembrane.log_stderr("# TMBETA-NET: Predicting strands for %s - %s\n" \
                      % (seqid, proteins[seqid]['name']))
    out = show()
    time.sleep(1)
    
    # parse the web page returned, extract strand boundaries
    proteins[seqid]['tmbeta_strands'] = []
    for l in out.split('\n'):
      if inmembrane.LOG_DEBUG: inmembrane.log_stderr("## " + l)

      if "<BR>Segment " in l:
        i,j = l.split(":")[1].split("to")
        i = int(i.strip()[1:])
        j = int(j.strip()[1:])
        proteins[seqid]['tmbeta_strands'].append([i,j])

        if inmembrane.LOG_DEBUG: inmembrane.log_stderr("# TMBETA-NET(web) segments: %s, %s" % (i, j))

    tmbeta_strands[seqid] = proteins[seqid]['tmbeta_strands']

  # we store the parsed strand boundaries in JSON format
  fh = open(outfile, 'w')
  fh.write(json.dumps(tmbeta_strands, separators=(',',':\n')))
  fh.close()

  return tmbeta_strands
