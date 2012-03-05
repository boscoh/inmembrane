import os, time, StringIO

import twill
from twill.commands import find, formfile, follow, fv, go, show, \
                             showforms, showlinks, submit
                             
from inmembrane import error_output, parse_fasta_header

__DEBUG__ = False

def bomp_web(params, proteins, \
             url="http://services.cbu.uib.no/tools/bomp/", force=False):
  """
  Uses the BOMP web service (http://services.cbu.uib.no/tools/bomp/) to
  predict if proteins are outer membrane beta-barrels.
  """
  
  bomp_out = 'bomp.out'
  error_output("# BOMP(web) %s > %s" % (params['fasta'], bomp_out))
  
  if not force and os.path.isfile(bomp_out):
    error_output("# -> skipped: %s already exists" % bomp_out)
    bomp_categories = {}
    fh = open(bomp_out, 'r')
    for l in fh:
      words = l.split()
      bomp_category = int(words[-1:][0])
      seqid = parse_fasta_header(l)[0]
      proteins[seqid]['bomp'] = bomp_category
      bomp_categories[seqid] = bomp_category
    fh.close()
    return bomp_categories
  
  # dump extraneous output into this blackhole so we don't see it
  if not __DEBUG__: twill.set_output(StringIO.StringIO())
  
  go(url)
  if __DEBUG__: showforms()
  formfile("1", "queryfile", params["fasta"])
  submit()
  if __DEBUG__: show()
  
  # extract the job id from the page
  links = showlinks()
  job_id = None
  for l in links:
    if l.url.find("viewOutput") != -1:
      # grab job id from "viewOutput?id=16745338"
      job_id = int(l.url.split("=")[1])
  
  if __DEBUG__: print "BOMP job id: ", job_id
  
  if not job_id:
    # something went wrong
    error_output("# BOMP error: Can't find job id")
    return
  
  # parse the HTML table and extract categories
  go("viewOutput?id=%i" % (job_id))
  
  polltime = 10
  error_output("# Waiting for BOMP to finish .")
  while True:
    try:
      find("Not finished")
      error_output(".")
    except:
      # Finished ! Pull down the result page.
      error_output(". done!\n")
      go("viewOutput?id=%i" % (job_id))
      if __DEBUG__: print show()
      break
      
    # Not finished. We keep polling for a time until
    # we give up
    time.sleep(polltime)
    polltime = polltime * 2
    if polltime >= 7200: # 2 hours
      error_output("# BOMP error: Taking too long.")
      return
    go("viewOutput?id=%i" % (job_id))
    if __DEBUG__: print show()
      
  bomp_html = show()
  if __DEBUG__: print bomp_html
  
  # Results are in the only <table> on this page, formatted like:
  # <tr><th>gi|107836852|gb|ABF84721.1<th>5</tr>
  from BeautifulSoup import BeautifulSoup
  soup = BeautifulSoup(bomp_html)
  bomp_categories = {} # dictionary of {name, category} pairs
  for tr in soup.findAll('tr')[1:]:
    n, c = tr.findAll('th')
    name = parse_fasta_header(n.text.strip())[0]
    category = int(c.text)
    bomp_categories[name] = category
  
  # write BOMP results to a tab delimited file
  fh = open(bomp_out, 'w')
  for k,v in bomp_categories.iteritems():
    fh.write("%s\t%i\n" % (k,v))
  fh.close()
  
  if __DEBUG__: print bomp_categories
  
  # label proteins with bomp classification (int) or False
  for name in proteins:
    if "bomp" not in proteins[name]:
      if name in bomp_categories:
        category = int(bomp_categories[name])
        proteins[name]['bomp'] = category
      else:
        proteins[name]['bomp'] = False
  
  if __DEBUG__: print proteins
  
  return bomp_categories
  
  """
  # Alternative: just get binary classification results via the
  #              FASTA output BOMP links to
  #
  # use the job id to jump straight to the fasta results
  # if a sequence is here, it's classified as an OMP barrel
  go("viewFasta?id=%i" % (job_id))
  bomp_seqs = show()
  bomp_fasta_headers = read_fasta_keys(StringIO.StringIO(show()))
  # label the predicted TMBs
  for name in bomp_fasta_headers:
    proteins[name]['bomp'] = True
    
  # label all the non-TMBs
  for name in proteins:
    if "bomp" not in proteins[name]:
      proteins[name]['bomp'] = False
  """
