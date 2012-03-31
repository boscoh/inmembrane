import sys, os, time, StringIO

import twill
from twill.commands import find, formfile, follow, fv, go, show, \
                             showforms, showlinks, submit, agent
import helpers

def parse_tatfind_output(output, proteins):
  """
  Parses the TatFind HTML output (file-like object or a list of strings) 
  an uses it to annotate and return an associated 'proteins' data structure.
  """
  for l in output:
    if "Results for" in l:
      seqid = l.split("Results for ")[1].split(":")[:-1][0]
      # parse id string to bring it to our format
      seqid, unused = helpers.parse_fasta_header(seqid)
      # "TRUE" or "FALSE"
      tat_pred = l.split("Results for ")[1].split(":")[-1:][0].strip()
      if tat_pred == "TRUE":
        proteins[seqid]["is_tatfind"] = True
      else:
        proteins[seqid]["is_tatfind"] = False
        
  return proteins

def annotate_tatfind_web(params, proteins, \
             url="http://signalfind.org/tatfind.html", force=False):
  """
  Interfaces with the TatFind web service at (http://signalfind.org/tatfind.html) 
  to predict if protein sequences contain Twin-Arginine Translocation (Tat) 
  signal peptides.
  """
  # set the user-agent so web services can block us if they want ... :/
  python_version = sys.version.split()[0]
  agent("Python-urllib/%s (twill; inmembrane)" % python_version)
  
  outfn = 'tatfind.out'
  helpers.log_stderr("# TatFind(web) %s > %s" % (params['fasta'], outfn))
  
  if not force and os.path.isfile(outfn):
    helpers.log_stderr("# -> skipped: %s already exists" % outfn)
    fh = open(outfn, 'r')
    proteins = parse_tatfind_output(fh, proteins)
    fh.close()
    return proteins
  
  # dump extraneous output into this blackhole so we don't see it
  if not helpers.LOG_DEBUG: twill.set_output(StringIO.StringIO())
  
  go(url)
  if helpers.LOG_DEBUG: showforms()
  formfile("1", "seqFile", params["fasta"])
  submit()
  if helpers.LOG_DEBUG: show()
  
  tatfind_output = show()
  if helpers.LOG_DEBUG: helpers.log_stderr(tatfind_output)
  
  # write raw TatFind output to a file
  fh = open(outfn, 'w')
  fh.write(tatfind_output)
  fh.close()
  
  proteins = parse_tatfind_output(tatfind_output.split("\n"), proteins) 
  
  return proteins
