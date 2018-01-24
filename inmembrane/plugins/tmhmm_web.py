# -*- coding: utf-8 -*-

citation = {'ref': u"Anders Krogh, Björn Larsson, Gunnar von Heijne and Erik "
                   u"L. L. Sonnhammer (2001) Predicting Transmembrane Protein "
                   u"Topology with a Hidden Markov Model: Application to "
                   u"Complete Genomes. J. Mol. Biol. 305:567-580. \n"
                   u"<http://dx.doi.org/10.1006/jmbi.2000.4315>",
            'name': "TMHMM 2.0"
            }

__DEBUG__ = False

import sys, os, time, json, re
from suds.client import Client
from suds.bindings import binding
import logging
from inmembrane.helpers import log_stderr, print_proteins


def annotate(params, proteins, \
             # url = 'http://www.cbs.dtu.dk/ws/TMHMM/TMHMM_2_0_ws0.wsdl', \
             # url = 'http://www.cbs.dtu.dk/ws/TMHMM/TMHMM_2_0_ws1.wsdl', \
             url="http://raw.github.com/boscoh/inmembrane/master/inmembrane/plugins/extra/TMHMM_2_0_ws0.wsdl", \
             batchsize=2000, \
             force=False):
    mapping = {'TMhelix': 'tmhmm_helices', \
               'outside': 'tmhmm_outer_loops', \
               'inside': 'tmhmm_inner_loops', \
               }

    if __DEBUG__:
        logging.basicConfig(level=logging.INFO)
        # soap messages (in&out) and http headers
        logging.getLogger('suds.client').setLevel(logging.DEBUG)

        # grab the cached results if present
    outfile = "tmhmm_web.out"
    if not force and os.path.isfile(outfile):
        log_stderr("# -> skipped: %s already exists" % outfile)
        fh = open(outfile, 'r')
        annots = json.loads(fh.read())
        fh.close()
        for seqid in annots:
            for k in mapping.values():
                proteins[seqid][k] = annots[seqid][k]

        citation['name'] = annots[seqid]['program_name']
        return proteins

    log_stderr("# TMHMM(web), %s > %s" % (params['fasta'], outfile))
    log_stderr(
        "# TMHMM(web): submitting in batches of %i sequences" % batchsize)

    seqids = proteins.keys()
    tmhmm_dict = {}
    while seqids:
        seqid_batch = seqids[0:batchsize]
        del seqids[0:batchsize]
        client = Client(url, cache=None)
        request = client.factory.create('runService.parameters')

        # this is a horrible horrible workaround to account for the fact that
        # the lipop SOAP service returns null results if there is are certain
        # non-alphanumeric characters in the sequence id provided. horrible.
        tmhmm_seq_id_mapping = {}
        seqcount = 0

        sys.stderr.write("# ")
        for seqid in seqid_batch:
            seq = client.factory.create(
                'runService.parameters.sequencedata.sequence')

            # workaround: removes any non-alphanumeric character (except '_') and adds
            # a unique number to the start to ensure every id is unique after mangling
            newseqid = `seqcount` + re.sub(r'[^\w]', "", seqid)
            seqcount += 1
            tmhmm_seq_id_mapping[newseqid] = seqid
            # seq.id = seqid
            seq.id = newseqid

            seq.seq = proteins[seqid]['seq']
            request.sequencedata.sequence.append(seq)
            sys.stderr.write(".")

        response = client.service.runService(request)

        sys.stderr.write("\n")

        # pollQueue
        job = client.factory.create('pollQueue.job')
        job.jobid = response.jobid
        response = client.service.pollQueue(job)
        retries = 0
        sys.stderr.write("# Waiting for TMHMM(web) results ")
        while response.status != "FINISHED" and retries < 100:
            response = client.service.pollQueue(job)
            time.sleep(10 + (retries * 2))
            retries += 1
            sys.stderr.write(".")

            # if something goes wrong, note it and skip TMHMM
            # by returning
            if response.status == "REJECTED" or \
                    response.status == "UNKNOWN JOBID" or \
                    response.status == "QUEUE DOWN" or \
                    response.status == "FAILED":
                log_stderr("TMHMM(web) failed: '%s'" % (response.status))
                return proteins

        sys.stderr.write(" done !\n")

        # fetchResults
        done_job = client.factory.create('fetchResult.job')
        done_job.jobid = response.jobid
        result = client.service.fetchResult(done_job)

        if __DEBUG__: log_stderr(str(result))

        citation["name"] = result[0].method + " " + result[0].version

        for res in result.ann:
            # seqid = res.sequence.id
            seqid = tmhmm_seq_id_mapping[res.sequence.id]
            if 'tmhmm_helices' not in proteins[seqid]:
                proteins[seqid].update({
                    'tmhmm_helices': [],
                    'tmhmm_inner_loops': [],
                    'tmhmm_outer_loops': []
                })
            if len(res.annrecords) > 0:
                for segment in res.annrecords.annrecord:
                    if segment.comment in mapping:
                        tmhmmkey = mapping[segment.comment]
                        proteins[seqid][tmhmmkey].append( \
                            (segment.range.begin, segment.range.end))

            # for caching in the outfile
            if seqid not in tmhmm_dict:
                tmhmm_dict[seqid] = {}

            # extract a copy of results from proteins dictionary
            # ready to we written to cache file
            for k in mapping.values():
                tmhmm_dict[seqid][k] = proteins[seqid][k]

            tmhmm_dict[seqid]['program_name'] = citation['name']

    if __DEBUG__: print_proteins(proteins)
    # we store the minimal stuff in JSON format
    fh = open(outfile, 'w')
    fh.write(json.dumps(tmhmm_dict, separators=(',', ':\n')))
    fh.close()

    return proteins
