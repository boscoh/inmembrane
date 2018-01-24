# -*- coding: utf-8 -*-

citation = {'ref': u"Garrow, A.G., Agnew, A. and Westhead, D.R. TMB-Hunt: An "
                   u"amino acid composition based method to screen proteomes "
                   u"for beta-barrel transmembrane proteins. BMC "
                   u"Bioinformatics, 2005, 6: 56 \n "
                   u"<http://dx.doi.org/10.1186/1471-2105-6-56>",
            'name': "TMB-HUNT"
            }

__DEBUG__ = False

import sys, os, time, StringIO

import twill
from twill.commands import find, formfile, follow, fv, go, show, \
    showforms, showlinks, submit, agent

from inmembrane.helpers import log_stderr, parse_fasta_header


def annotate(params, proteins, \
             force=False):
    """
    DEPRECATED: The TMB-HUNT server appears to be permanently offline.

    Uses the TMB-HUNT web service
    (http://bmbpcu36.leeds.ac.uk/~andy/betaBarrel/AACompPred/aaTMB_Hunt.cgi) to
    predict if proteins are outer membrane beta-barrels.

    NOTE: In my limited testing, TMB-HUNT tends to perform very poorly in
          terms of false positives and false negetives. I'd suggest using only
          BOMP.
    """
    # TODO: automatically split large sets into multiple jobs
    #       TMB-HUNT will only take 10000 seqs at a time
    if len(proteins) >= 10000:
        log_stderr(
            "# ERROR: TMB-HUNT(web): can't take more than 10,000 sequences.")
        return

    # set the user-agent so web services can block us if they want ... :/
    python_version = sys.version.split()[0]
    agent("Python-urllib/%s (twill; inmembrane)" % python_version)

    out = 'tmbhunt.out'
    log_stderr("# TMB-HUNT(web) %s > %s" % (params['fasta'], out))

    if not force and os.path.isfile(out):
        log_stderr("# -> skipped: %s already exists" % out)
        return parse_tmbhunt(proteins, out)

    # dump extraneous output into this blackhole so we don't see it
    if not __DEBUG__: twill.set_output(StringIO.StringIO())

    go("http://bmbpcu36.leeds.ac.uk/~andy/betaBarrel/AACompPred/aaTMB_Hunt.cgi")
    if __DEBUG__: showforms()

    # read up the FASTA format seqs
    fh = open(params['fasta'], 'r')
    fasta_seqs = fh.read()
    fh.close()

    # fill out the form
    fv("1", "sequences", fasta_seqs)

    submit()
    if __DEBUG__: showlinks()

    # small jobs will lead us straight to the results, big jobs
    # go via a 'waiting' page which we skip past if we get it
    job_id = None
    try:
        # we see this with big jobs
        result_table_url = follow(
            "http://www.bioinformatics.leeds.ac.uk/~andy/betaBarrel/AACompPred/tmp/tmp_output.*.html")
        job_id = result_table_url.split('tmp_output')[-1:][0].split('.')[0]
    except:
        # small jobs take us straight to the html results table
        pass

    # parse the job_id from the url, since due to a bug in
    # TMB-HUNT the link on the results page from large jobs is wrong
    if not job_id: job_id = \
    follow("Full results").split('/')[-1:][0].split('.')[0]
    log_stderr(
        "# TMB-HUNT(web) job_id is: %s <http://www.bioinformatics.leeds.ac.uk/~andy/betaBarrel/AACompPred/tmp/tmp_output%s.html>" % (
        job_id, job_id))

    # polling until TMB-HUNT finishes
    # TMB-HUNT advises that 4000 sequences take ~10 mins
    # we poll a little faster than that
    polltime = (len(proteins) * 0.1) + 2
    while True:
        log_stderr("# TMB-HUNT(web): waiting another %i sec ..." % (polltime))
        time.sleep(polltime)
        try:
            go(
                "http://bmbpcu36.leeds.ac.uk/~andy/betaBarrel/AACompPred/tmp/%s.txt" % (
                job_id))
            break
        except:
            polltime = polltime * 2

        if polltime >= 7200:  # 2 hours
            log_stderr("# TMB-HUNT error: Taking too long.")
            return

    txt_out = show()

    # write raw TMB-HUNT results
    fh = open(out, 'w')
    fh.write(txt_out)
    fh.close()

    return parse_tmbhunt(proteins, out)


def parse_tmbhunt(proteins, out):
    """
    Takes the filename of a TMB-HUNT output file (text format)
    & parses the outer membrane beta-barrel predictions into the proteins dictionary.
    """
    # parse TMB-HUNT text output
    tmbhunt_classes = {}
    for l in open(out, 'r'):
        # inmembrane.log_stderr("# TMB-HUNT raw: " + l[:-1])
        if l[0] == ">":
            # TMB-HUNT munges FASTA ids by making them all uppercase,
            # so we find the equivalent any-case id in our proteins list
            # and use that. ugly but necessary.
            seqid, desc = parse_fasta_header(l)
            for i in proteins.keys():
                if seqid.upper() == i.upper():
                    seqid = i
                    desc = proteins[i]['name']

            probability = None
            classication = None
            tmbhunt_classes[seqid] = {}
        if l.find(
                "Probability of a NON-BETA BARREL protein with this score:") != -1:
            # we convert from probability of NON-BARREL to probability of BARREL
            probability = 1 - float(l.split(":")[1].strip())
        if l[0:11] == "Conclusion:":
            classication = l.split(":")[1].strip()
            if classication == "BBMP":
                tmbhunt_classes[seqid]['tmbhunt'] = True
                tmbhunt_classes[seqid]['tmbhunt_prob'] = probability

                proteins[seqid]['tmbhunt'] = True
                proteins[seqid]['tmbhunt_prob'] = probability

            elif classication == "Non BBMP":
                tmbhunt_classes[seqid]['tmbhunt'] = False
                tmbhunt_classes[seqid]['tmbhunt_prob'] = probability

                proteins[seqid]['tmbhunt'] = False
                proteins[seqid]['tmbhunt_prob'] = probability

    # inmembrane.log_stderr(str(tmbhunt_classes))
    return tmbhunt_classes
