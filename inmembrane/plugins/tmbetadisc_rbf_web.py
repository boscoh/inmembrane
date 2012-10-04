# -*- coding: utf-8 -*-

citation = {'ref': "Ou Y-YY, Gromiha MMM, Chen S-AA, Suwa M (2008) "
                   "TMBETADISC-RBF: Discrimination of beta-barrel "
                   "membrane proteins using RBF networks and PSSM profiles. "
                   "Computational biology and chemistry. \n"
                   "<http://dx.doi.org/10.1016/j.compbiolchem.2008.03.002>",
            'name': "TMBETADISC-RBF"
           }

__DEBUG__ = False

import sys, os, time, StringIO

import twill
from twill.commands import find, formfile, follow, fv, go, show, \
                             showforms, showlinks, submit, agent
from BeautifulSoup import BeautifulSoup                             
from inmembrane.helpers import log_stderr, parse_fasta_header, dict_get

def parse_tmbetadisc_output(output, proteins):
  """
  Parses the TMBETADISC-RBF output (file-like object or a list of strings) 
  an uses it to annotate and return an associated 'proteins' data structure.
  """

  soup = BeautifulSoup(output)
  # parse the table. we pop of single data cells one at a time
  fields = soup.findAll("td")
  fields.reverse()
  f = fields.pop() # discard first <td>1</td> field
  try:
    while len(fields) > 0:
      f = fields.pop().text
      seqid, result = parse_fasta_header(f)
      if "Non-Outer Membrane Protein" in result:
        proteins[seqid]["is_tmbetadisc_rbf"] = False
      elif "is Outer Membrane Protein" in result:
        proteins[seqid]["is_tmbetadisc_rbf"] = True
      fields.pop()
  except IndexError:
    # we get here when we run out of table fields to pop
    pass
        
  return proteins

def annotate(params, proteins, \
             url="http://rbf.bioinfo.tw/"+
                 "~sachen/OMPpredict/"+
                 "TMBETADISC-RBF-Content.html", force=False):
  """
  Interfaces with the TatFind web service at 
  (http://rbf.bioinfo.tw/~sachen/OMPpredict/TMBETADISC-RBF.php) 
  to predict if protein sequence is likely to be an outer membrane beta-barrel.
  
  Note that the default URL we use it different to the regular form used
  by web browsers, since we need to bypass some AJAX fun.
  """
  # TODO: automatically split large sets into multiple jobs
  #       since TMBETADISC seems to not like more than take 
  #       ~5000 seqs at a time
  if len(proteins) >= 5000:
    log_stderr("# ERROR: TMBETADISC-RBF(web): tends to fail with > ~5000 sequences.")
    return
  
  # set the user-agent so web services can block us if they want ... :/
  python_version = sys.version.split()[0]
  agent("Python-urllib/%s (twill; inmembrane)" % python_version)
  
  outfn = 'tmbetadisc-rbf.out'
  log_stderr("# TMBETADISC-RBF(web) %s > %s" % (params['fasta'], outfn))
  
  if not force and os.path.isfile(outfn):
    log_stderr("# -> skipped: %s already exists" % outfn)
    fh = open(outfn, 'r')
    proteins = parse_tmbetadisc_output(fh.read(), proteins)
    fh.close()
    return proteins
  
  # dump extraneous output into this blackhole so we don't see it
  if not __DEBUG__: twill.set_output(StringIO.StringIO())
  
  go(url)
  if __DEBUG__: showforms()
  formfile("1", "userfile", params["fasta"])
  fv("1", "format", "file")

  # set the user defined method
  method_map = {"aa":"Amino Acid Composition",
                "dp":"Depipetide Composition",
                "aadp":"Amino Acid & Depipetide Composition",
                "pssm":"PSSM"}
  if dict_get(params, 'tmbetadisc_rbf_method'):
    try:
      method = method_map[params['tmbetadisc_rbf_method']]
    except KeyError:
      log_stderr("# ERROR: Invalid setting from tmbetadisc_rbf_method. \
                    Must be set to aa, dp, aadp or pssm.")
      sys.exit()

  #fv("1", "select", "Amino Acid Composition")
  #fv("1", "select", "Depipetide Composition")
  #fv("1", "select", "Amino Acid & Depipetide Composition")
  #fv("1", "select", "PSSM")
  fv("1", "select", method)
  
  submit()
  
  waiting_page = show()
  if __DEBUG__: log_stderr(waiting_page)

  for l in waiting_page.split('\n'):
    if l.find("TMBETADISC-RBF-action.php?UniqueName=") != -1:
      result_url = l.split("'")[1]

  time.sleep(5)
  
  go(result_url)
  
  output = show()
  if __DEBUG__: log_stderr(output)
  
  # write raw output to a file
  fh = open(outfn, 'w')
  fh.write(output)
  fh.close()
  
  proteins = parse_tmbetadisc_output(output, proteins) 
  
  return proteins
