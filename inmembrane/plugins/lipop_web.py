# -*- coding: utf-8 -*-
# code loosely based on: https://portal.nordu.net/display/ndgfwiki/Python-Suds
# see: http://www.cbs.dtu.dk/ws/ws.php?entry=SignalP4

citation = {'ref': u"Agnieszka S. Juncker, Hanni Willenbrock, "
                   u"Gunnar Von Heijne, Søren Brunak, Henrik Nielsen, "
                   u"And Anders Krogh. (2003) Prediction of lipoprotein "
                   u"signal peptides in Gram-negative bacteria. Protein "
                   u"Science 12:1652–1662. \n"
                   u"<http://dx.doi.org/10.1110/ps.0303703>",
            'name': "LipoP 1.0.ws0"
            }

__DEBUG__ = False

import sys, os, time, json, re, urllib2
from suds.client import Client
from suds.bindings import binding
import logging
from inmembrane.helpers import log_stderr


def annotate(params, proteins, \
             # url = 'http://www.cbs.dtu.dk/ws/LipoP/LipoP_1_0_ws0.wsdl', \
             # we host our own fixed version of the WSDL for the moment
             url="http://raw.github.com/boscoh/inmembrane/master/inmembrane/plugins/extra/LipoP_1_0_ws0.wsdl", \
             # url = "http://www.cbs.dtu.dk/ws/LipoP/LipoP_1_0_ws0.wsdl",
             batchsize=2000, \
             force=False):
    if __DEBUG__:
        logging.basicConfig(level=logging.INFO)
        # soap messages (in&out) and http headers
        logging.getLogger('suds.client').setLevel(logging.DEBUG)

        # grab the cached results if present
    outfile = "lipop_web.out"
    if not force and os.path.isfile(outfile):
        log_stderr("# -> skipped: %s already exists" % outfile)
        fh = open(outfile, 'r')
        annots = json.loads(fh.read())
        fh.close()
        for seqid in annots:
            proteins[seqid]['is_lipop'] = annots[seqid]['is_lipop']
            proteins[seqid]['lipop_cleave_position'] = \
                annots[seqid]['lipop_cleave_position']

        citation['name'] = annots[seqid]['program_name']
        return proteins

    log_stderr("# LipoP(web), %s > %s" % (params['fasta'], outfile))
    log_stderr(
        "# LipoP(web): submitting in batches of %i sequences" % batchsize)

    """
    # ensure schemas are correctly imported (workaround for broken schemas ..)
    from suds.xsd.doctor import ImportDoctor
    from suds.xsd.doctor import Import
    imp = Import("http://www.cbs.dtu.dk/ws/ws-common", location="http://www.cbs.dtu.dk/ws/common/ws_common_1_0b.xsd")
    imp.filter.add("http://www.cbs.dtu.dk/ws/WSLipoP_1_0_ws0")
    doctor = ImportDoctor(imp)
    client = Client(url, doctor=doctor, cache=None)
    #client = Client(url, plugins=[doctor], cache=None)
    """

    seqids = proteins.keys()
    lipop_dict = {}
    while seqids:
        seqid_batch = seqids[0:batchsize]
        del seqids[0:batchsize]

        client = Client(url, cache=None)

        request = client.factory.create('runService.parameters')

        # this is a horrible horrible workaround to account for the fact that
        # the lipop SOAP service returns null results if there is are certain
        # non-alphanumeric characters in the sequence id provided. horrible.
        lipop_seq_id_mapping = {}
        seqcount = 0

        sys.stderr.write("# ")
        for seqid in seqid_batch:
            seq = client.factory.create(
                'runService.parameters.sequencedata.sequence')

            # workaround: removes any non-alphanumeric character (except '_') and adds
            # a unique number to the start to ensure every id is unique after mangling
            newseqid = `seqcount` + re.sub(r'[^\w]', "", seqid)
            seqcount += 1
            lipop_seq_id_mapping[newseqid] = seqid
            seq.id = newseqid
            # seq.id = seqid
            seq.seq = proteins[seqid]['seq']

            request.sequencedata.sequence.append(seq)
            sys.stderr.write(".")
        try:
            response = client.service.runService(request)
        except urllib2.URLError as e:
            log_stderr("ERROR LipoP(web) failed: '%s'" % `e.reason`)
            return proteins

        sys.stderr.write("\n")

        # pollQueue
        job = client.factory.create('pollQueue.job')
        job.jobid = response.jobid
        response = client.service.pollQueue(job)
        retries = 0
        sys.stderr.write("# Waiting for LipoP(web) results ")
        while response.status != "FINISHED" and retries < 12:
            response = client.service.pollQueue(job)
            time.sleep(10 + (retries ** 2))
            retries += 1
            sys.stderr.write(".")

            # if something goes wrong, note it and skip LipoP
            # by returning
            if response.status == "REJECTED" or \
                    response.status == "UNKNOWN JOBID" or \
                    response.status == "QUEUE DOWN" or \
                    response.status == "FAILED":
                log_stderr("LipoP(web) failed: '%s'" % (response.status))
                return proteins

        sys.stderr.write(" done !\n")

        # fetchResults
        done_job = client.factory.create('fetchResult.job')
        done_job.jobid = response.jobid
        result = client.service.fetchResult(done_job)
        if __DEBUG__: log_stderr(str(result))

        citation["name"] = result[0].method + " " + result[0].version

        # TODO: the better way to do this would be to save the entire SOAP
        #       response returned by client.last_received() and then parse
        #       that upon plugin invocation (above) using suds.sax
        #       This way we save everything in the analysis, not just
        #       the details we are interested in right now
        for res in result.ann:
            # seqid = res.sequence.id
            seqid = lipop_seq_id_mapping[res.sequence.id]
            # init as if no lipop hit, may be reset below
            proteins[seqid]['is_lipop'] = False
            proteins[seqid]['lipop_cleave_position'] = 0
            proteins[seqid]['lipop_im_retention_signal'] = False
            if len(res.annrecords) > 0:
                # range.end - this is the first residue (Cys) of the mature protein if
                #  there is a SpII cleavage site
                for annrec in res.annrecords.annrecord:
                    if annrec.feature == "CleavII":
                        proteins[seqid]['lipop_cleave_position'] = int(
                            annrec.range.begin)
                        proteins[seqid]['is_lipop'] = True

                        # check for an E.coli style inner membrane retention signal
                        # Asp+2 to cleavage site. There are other apparent retention
                        # signals in E. coli and other gram- bacteria in addition to
                        # the Asp+2 which we don't detect here (yet).
                        # (Yamaguchi et al, 1988; Tokuda and Matsuyama, 2005 [review])
                        plus_two = proteins[seqid]['lipop_cleave_position'] + 1
                        if proteins[seqid]['seq'][plus_two] == 'D':
                            proteins[seqid]['lipop_im_retention_signal'] = True

            # for caching in the outfile
            if seqid not in lipop_dict:
                lipop_dict[seqid] = {}
            lipop_dict[seqid]['is_lipop'] = proteins[seqid]['is_lipop']
            lipop_dict[seqid]['lipop_cleave_position'] = \
                proteins[seqid]['lipop_cleave_position']
            lipop_dict[seqid]['lipop_im_retention_signal'] = \
                proteins[seqid]['lipop_im_retention_signal']
            lipop_dict[seqid]['program_name'] = citation['name']

    # we store the minimal stuff in JSON format
    fh = open(outfile, 'w')
    fh.write(json.dumps(lipop_dict, separators=(',', ':\n')))
    fh.close()

    return proteins
