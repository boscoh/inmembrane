# code loosely based on: https://portal.nordu.net/display/ndgfwiki/Python-Suds
# see: http://www.cbs.dtu.dk/ws/ws.php?entry=SignalP4

citation = {'ref': u"Petersen TN, Brunak S, von Heijne G, Nielsen H. "
                   u"SignalP 4.0: discriminating signal peptides from "
                   u"transmembrane regions. Nature methods 2011 "
                   u"Jan;8(10):785-6. \n"
                   u"<http://dx.doi.org/10.1038/nmeth.1701>",
            'name': "SignalP 4.0"
           }
          
__DEBUG__ = False
import sys, os, time, json
from suds.client import Client
from suds.bindings import binding 
import logging
from helpers import log_stderr

if __DEBUG__:
  logging.basicConfig(level=logging.INFO)
  # soap messages (in&out) and http headers
  logging.getLogger('suds.client').setLevel(logging.DEBUG)

def annotate(params, proteins, \
             #url = 'http://www.cbs.dtu.dk/ws/SignalP4/SignalP4_4_0_ws0.wsdl', \
             url = 'http://www.cbs.dtu.dk/ws/SignalP/SignalP_3_1_ws0.wsdl', \
             force=False):
             
  
  # TODO: automatically split large sets into multiple jobs
  #       since the SignalP webservice will take a maximum of
  #       2000 sequences at a time
  if len(proteins) > 2000:
    log_stderr("# ERROR: SignalP(web): currently only works with < 2000 sequences.")
    return
    
  # grab the cached results of present
  outfile = "signalp_web.out"
  if not force and os.path.isfile(outfile):
    log_stderr("# -> skipped: %s already exists" % out)
    fh = open(outfile, 'r')
    annot = json.loads(fh.read())
    fh.close()
    for seqid in annots:
      proteins[seqid]['is_signalp'] = annots[seqid]['is_signalp']
      proteins[seqid]['signalp_cleave_position'] = \
        annots[seqid]['signalp_cleave_position']
        
      return proteins
  
  log_stderr("# SignalP(web), %s > %s" % (params['fasta'], outfile))
  
  client = Client(url, cache=None)
  request=client.factory.create('runService.parameters')
  
  sys.stderr.write("# ")
  for seqid in proteins:
    seq = client.factory.create('runService.parameters.sequencedata.sequence')
    seq.id = seqid
    seq.seq = proteins[seqid]['seq']
  
    # organism can be 'euk', 'gram+', 'siganlgram-'
    request.organism = params['signalp4_organism']
    # default for SignalP 4.0
    #request.method = 'best'
    # default for SignalP 3.1
    #request.method = 'nn+hmm'
    request.sequencedata.sequence.append(seq)
  
    response = client.service.runService(request)
    sys.stderr.write(".")
  sys.stderr.write("\n")
  
  #pollQueue
  job = client.factory.create('pollQueue.job')
  job.jobid = response.jobid
  response = client.service.pollQueue(job)
  retries = 0
  sys.stderr.write("# Waiting for SignalP results ")
  while response.status != "FINISHED" and retries < 100:
    response = client.service.pollQueue(job)
    time.sleep(10 + (retries*2))
    retries += 1
    sys.stderr.write(".")
    
    # if something goes wrong, note it and skip SignalP
    # by returning
    if response.status == "REJECTED" or \
       response.status == "UNKNOWN JOBID" or \
       response.status == "QUEUE DOWN":
      log_stderr("SignalP(web) failed: '%s'" % (response.status))
      return proteins
      
  sys.stderr.write(" done !\n")
    
  #fetchResults
  done_job = client.factory.create('fetchResult.job')
  done_job.jobid = response.jobid
  result = client.service.fetchResult(done_job)
  
  # end of signal-nn
  
  citation["name"] = result[0].method + " " + result[0].version
    
  # TODO: the better way to do this would be to save the entire SOAP
  #       response returned by client.last_received() and then parse
  #       that upon plugin invocation (above) using suds.sax
  #       This way we save everything in the analysis, not just
  #       the details we are interested in right now
  signalp_dict = {}
  for res in result.ann:
    seqid = res.sequence.id
    # range.end - this is the first residue of the mature protein if
    #  there is a cleavage site
    cleavage_site = int(res.annrecords.annrecord[0].range.end)
    proteins[seqid]['signalp_cleave_position'] = cleavage_site
    # from 'comment', "Y" or "N noTm" or "N TM" where "Y" means signal peptide
    signal_yn = res.annrecords[0][0].comment[0]
    if signal_yn == "Y":
      proteins[seqid]['is_signalp'] = True
    else:
      proteins[seqid]['is_signalp'] = False
      
    # for caching in the outfile
    if seqid not in signalp_dict:
      signalp_dict[seqid] = {}
    signalp_dict[seqid]['is_signalp'] = proteins[seqid]['is_signalp']
    signalp_dict[seqid]['signalp_cleave_position'] = \
      proteins[seqid]['signalp_cleave_position']
  
  # we store the minimal stuff in JSON format
  fh = open(outfile, 'w')
  fh.write(json.dumps(signalp_dict, separators=(',',':\n')))
  fh.close()    
  
  return proteins

