# code loosely based on: https://portal.nordu.net/display/ndgfwiki/Python-Suds
# see: http://www.cbs.dtu.dk/ws/ws.php?entry=TatP
__DEBUG__ = True
import sys, time
from suds.client import Client
from suds.bindings import binding 
import logging
from helpers import log_stderr
  
def annotate(params, proteins, \
             url="http://www.cbs.dtu.dk/ws/TatP/TatP_1_0_ws0.wsdl", force=False):
  if __DEBUG__:
    logging.basicConfig(level=logging.INFO)
    # SOAP messages (in&out) and http headers
    logging.getLogger('suds.client').setLevel(logging.DEBUG)
    
  outfile = "tatp_web.out"
  # TODO: check for cached results file
  
  log_stderr("# TatP(web), %s > %s" % (params['fasta'], outfile))
  
  # add schemas to work around for broken WSDL ...
  from suds.xsd.doctor import ImportDoctor
  from suds.xsd.doctor import Import
  imp = Import("http://www.cbs.dtu.dk/ws/ws-common", location='http://www.cbs.dtu.dk/ws/common/ws_common_1_0b.xsd')
  imp2 = Import("http://www.cbs.dtu.dk/ws/ws-tatp", location="http://www.cbs.dtu.dk/ws/TatP/ws_tatp_1_0_ws0.xsd")
  imp.filter.add("http://www.cbs.dtu.dk/ws/WSTatP_1_0_ws0")
  imp2.filter.add("http://www.cbs.dtu.dk/ws/WSTatP_1_0_ws0")
  doctor = ImportDoctor(imp, imp2)
  #client = Client(url, doctor=doctor, cache=None)
  client = Client(url, plugins=[doctor], cache=None)
  
  # Alternative suds example for adding schemas explicitly:
  #  http://stackoverflow.com/questions/4719854/soap-suds-and-the-dreaded-schema-type-not-found-error

  sys.stderr.write("# ")
  for seqid in proteins:
    seq = client.factory.create('runService.parameters.sequencedata.sequence')
    seq.id = seqid
    seq.seq = proteins[seqid]['seq']
  
    request=client.factory.create('runService.parameters')
    request.sequencedata.sequence=[seq]
    sys.stderr.write(".")
  
  response = client.service.runService(request)
  sys.stderr.write("\n")
  
  #pollQueue
  job = client.factory.create('pollQueue.job')
  job.jobid = response.jobid
  response = client.service.pollQueue(job)
  retries = 0
  sys.stderr.write("# Waiting for TatP results ")
  while response.status != "FINISHED" and retries < 100:
    response = client.service.pollQueue(job)
    time.sleep(10 + (retries*2))
    retries += 1
    sys.stderr.write(".")
    
    # if something goes wrong, note it and skip TatP
    # by returning
    if response.status == "REJECTED" or \
       response.status == "UNKNOWN JOBID" or \
       response.status == "QUEUE DOWN":
      log_stderr("TatP(web) failed: '%s'" % (response.status))
      return proteins
      
  sys.stderr.write(" done !\n")
  
  #fetchResults
  done_job = client.factory.create('fetchResult.job')
  done_job.jobid = response.jobid
  result = client.service.fetchResult(done_job)
  
  # TODO: this is returning as None.
  #       is the service actually broken ?
  print result
  print
  
