# -*- coding: utf-8 -*-
citation = {'ref': u"Agnieszka S. Juncker, Hanni Willenbrock, "
                   u"Gunnar Von Heijne, Søren Brunak, Henrik Nielsen, "
                   u"And Anders Krogh. (2003) Prediction of lipoprotein "
                   u"signal peptides in Gram-negative bacteria. Protein "
                   u"Science 12:1652–1662. \n"
                   u"<http://dx.doi.org/10.1110/ps.0303703>",
            'name': "LipoP 1.0 (web interface)"
           }
          
__DEBUG__ = False

import sys, os, time
from StringIO import StringIO
from BeautifulSoup import BeautifulSoup
import requests
from collections import OrderedDict
import inmembrane
from inmembrane.plugins.lipop1 import parse_lipop
from inmembrane.helpers import log_stderr
from inmembrane.helpers import generate_safe_seqids, proteins_to_fasta

def annotate(params, proteins, \
             batchsize = 2000, \
             force=False):
  baseurl = "http://www.cbs.dtu.dk"
  url = baseurl + "/cgi-bin/webface2.fcgi"

  if len(proteins) > 4000:
    log_stderr("# ERROR: The LipoP web inteface does not accept > \
                4000 sequences at one time.")
    return

  # grab the cached results if present
  outfile = "lipop_scrape_web.out"
  if not force and os.path.isfile(outfile):
    log_stderr("# -> skipped: %s already exists" % outfile)
    proteins, id_mapping = generate_safe_seqids(proteins)
    fh = open(outfile, 'r')
    resultpage = fh.read()
    fh.close()
    resultpage = resultpage.replace("<hr>", '\n')
    resultpage = resultpage.replace("<P>", '\n')
    soup = BeautifulSoup(resultpage)
    proteins = parse_lipop(soup('pre')[0].text, proteins, id_mapping=id_mapping)
    return proteins

  # get sequences in fasta format with munged ids (workaround for
  #  lipop sequence id munging)
  proteins, id_mapping = generate_safe_seqids(proteins)
  safe_fasta = proteins_to_fasta(proteins, use_safe_seqid=True)

  # we use an OrderedDict rather than a normal dictionary to work around 
  # some quirks in the CBS CGI (the server expects parameters in a certain 
  # order in the HTTP headers).
  payload = OrderedDict([('configfile',
                        "/usr/opt/www/pub/CBS/services/LipoP-1.0/LipoP.cf"),
                        ("SEQ",""),
                        ("outform","-noplot")])

  #files = {'seqfile': open(params['fasta'], 'rb')}
  files = {'seqfile': StringIO(safe_fasta)}

  log_stderr("# LipoP(scrape_web), %s > %s" % (params['fasta'], outfile))
  # log_stderr("# LipoP(scrape_web): submitting in batches of %i sequences" % batchsize)

  headers = {"User-Agent": 
             "python-requests/%s (inmembrane/%s)" % 
             (requests.__version__, inmembrane.__version__) }
  r = requests.post(url, data=payload, files=files, headers=headers)
  if __DEBUG__:
    print r.text
    # Example:
    #
    # <HTML>
    # <HEAD><TITLE>Webface Jobsubmission</TITLE></HEAD>
    # If Javascript is disabled, follow <a href="/cgi-bin/webface?jobid=LipoP,50B5432A10A9CD51&opt=wait">This link</a>
    #
    # <script LANGUAGE="JavaScript"><!--
    # location.replace("/cgi-bin/webface?jobid=LipoP,50B5432A10A9CD51&opt=wait")
    # //--></script>
    # </HTML>

  # extract the result URL (or die if job is rejected ...)
  if "Job rejected" in r.text:
    sys.stderr.write(r.text)
    sys.exit()
  soup = BeautifulSoup(r.text)
  resultlink = baseurl + soup.findAll('a')[0]['href']

  # brief pause (LipoP is quick), then grab the results at the result url
  sys.stderr.write("# Waiting briefly for LipoP results")
  time.sleep(len(proteins)/500)
  resultpage = requests.post(resultlink).text
  sys.stderr.write(" .. done !\n")

  if __DEBUG__:
    print resultpage
    # Example:
    #
    #   <pre>
    # # lcl_AE004092.1_cdsid_AAK33146.1 CYT score=-0.200913
    # # Cut-off=-3
    # lcl_AE004092.1_cdsid_AAK33146.1	LipoP1.0:Best	CYT	1	1	-0.200913
    # <P>
    # <hr>
    # # lcl_AE004092.1_cdsid_AAK33147.1 CYT score=-0.200913
    # # Cut-off=-3
    # lcl_AE004092.1_cdsid_AAK33147.1	LipoP1.0:Best	CYT	1	1	-0.200913
    # <P>
    # <hr>

  # replacing these tags make the output parsable by the
  # existing standalone lipop1 parser.
  resultpage = resultpage.replace("<hr>", '\n')
  resultpage = resultpage.replace("<P>", '\n')

  soup = BeautifulSoup(resultpage)
  proteins = parse_lipop(soup('pre')[0].text, proteins, id_mapping=id_mapping)

  # we store the downloaded page
  fh = open(outfile, 'w')
  fh.write(resultpage)
  fh.close()

  return proteins
