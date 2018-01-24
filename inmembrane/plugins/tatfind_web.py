# -*- coding: utf-8 -*-

citation = {'ref': u"Rose, R.W., T. Brüser,. J. C. Kissinger, and M. "
                   u"Pohlschröder. 2002. Adaptation of protein secretion "
                   u"to extremely high salt concentrations by extensive use "
                   u"of the twin arginine translocation pathway. Mol. Microbiol."
                   u"5: 943-950 \n"
                   u"<http://dx.doi.org/10.1046/j.1365-2958.2002.03090.x>",
            'name': "TatFind 1.4"
            }

__DEBUG__ = False

import sys, os, time, StringIO

import twill
from twill.commands import find, formfile, follow, fv, go, show, \
    showforms, showlinks, submit, agent
from inmembrane.helpers import log_stderr, parse_fasta_header


def parse_tatfind_output(output, proteins):
    """
    Parses the TatFind HTML output (file-like object or a list of strings)
    an uses it to annotate and return an associated 'proteins' data structure.
    """
    for l in output:
        if "Results for" in l:
            seqid = l.split("Results for ")[1].split(":")[:-1][0]
            # parse id string to bring it to our format
            seqid, unused = parse_fasta_header(seqid)
            # "TRUE" or "FALSE"
            tat_pred = l.split("Results for ")[1].split(":")[-1:][0].strip()
            if tat_pred == "TRUE":
                proteins[seqid]["is_tatfind"] = True
            else:
                proteins[seqid]["is_tatfind"] = False

    return proteins


def annotate(params, proteins, \
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
    log_stderr("# TatFind(web) %s > %s" % (params['fasta'], outfn))

    if not force and os.path.isfile(outfn):
        log_stderr("# -> skipped: %s already exists" % outfn)
        fh = open(outfn, 'r')
        proteins = parse_tatfind_output(fh, proteins)
        fh.close()
        return proteins

    # dump extraneous output into this blackhole so we don't see it
    if not __DEBUG__: twill.set_output(StringIO.StringIO())

    go(url)
    if __DEBUG__: showforms()
    formfile("1", "seqFile", params["fasta"])
    submit()
    if __DEBUG__: show()

    tatfind_output = show()
    if __DEBUG__: log_stderr(tatfind_output)

    # write raw TatFind output to a file
    fh = open(outfn, 'w')
    fh.write(tatfind_output)
    fh.close()

    proteins = parse_tatfind_output(tatfind_output.split("\n"), proteins)

    return proteins
