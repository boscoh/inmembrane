#
# Common helper functions for inmembrane
#

import inmembrane
import os, subprocess, sys, re
from collections import OrderedDict
import textwrap
from bs4 import BeautifulSoup

LOG_SILENT = False


def dict_get(this_dict, prop):
    """
    Useful helper function to get values from dicts. Takes in
    the possibility that the key may not exist, and in that
    case returns False. Makes it easier to write code to avoid
    handling the case of missing keys.
    """
    if prop not in this_dict:
        return False
    return this_dict[prop]


def run_with_output(cmd):
    """
    Runs an external program as a child process and captures
    the input. May not work with external programs that include
    complicated levels of shells.
    """
    p = subprocess.Popen(
        cmd, shell=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    return p.stdout.read()


def run(cmd, out_file=None):
    """
    Wrapper function to run external program so that existence
    of binary can be checked and the output directed to a specific
    file.
    """
    log_stderr("# " + cmd + " > " + out_file)
    if os.path.isfile(out_file) and (out_file != None):
        log_stderr("# -> skipped: %s already exists" % out_file)
        return
    binary = cmd.split()[0]
    is_binary_there = False
    if os.path.isfile(binary):
        is_binary_there = True
    if run_with_output('which ' + binary):
        is_binary_there = True
    if not is_binary_there:
        raise IOError("Couldn't find executable binary '" + binary + "'")
    if out_file:
        os.system(cmd + " > " + out_file)
    else:
        os.system(cmd)


def silence_log(b):
    """
    Turns on logging silent mode w.r.t. to log_stderr and log_stdout
    """
    global LOG_SILENT
    LOG_SILENT = b


def log_stderr(s, width=76, comment=True):
    """
    Wrapper for all stderr out. Allows future customization.
    """
    if LOG_SILENT:
        return
    if s and s[-1] != "\n":
        s += "\n"
    if not s.startswith("#"):
        s = "# " + s
    sys.stderr.write(s)


def log_stdout(s, width=76):
    """
    Wrapper for all stdout output. Allows future customization.
    """
    if LOG_SILENT:
        return
    print s


def parse_fasta_header(header):
    """
    Parses a FASTA format header (with our without the initial '>') and returns a
    tuple of sequence id and sequence name/description.

    If NCBI SeqID format (gi|gi-number|gb|accession etc, is detected
    the first id in the list is used as the canonical id (see see
    http://www.ncbi.nlm.nih.gov/books/NBK21097/#A631 ).
    """
    if header[0] == '>':
        header = header[1:]
    tokens = header.split('|')
    # check to see if we have an NCBI-style header
    if header.find("|") != -1 and len(tokens[0]) <= 3:
        # "gi|ginumber|gb|accession bla bla" becomes "gi|ginumber"
        seqid = "%s|%s" % (tokens[0], tokens[1].split()[0])
        name = tokens[-1:][0].strip()
    # otherwise just split on spaces & hope for the best
    else:
        tokens = header.split()
        seqid = tokens[0]
        name = header[0:-1].strip()

    return seqid, name


def seqid_to_filename(seqid):
    """
    Makes a sequence id filename friendly.
    (eg, replaces '|' with '_')
    """
    return seqid.replace("|", "_")


def create_proteins_dict(fasta):
    """
    The main data-structure used in inmembrane is the proteins dictionary,
    which is initialized here from a source FASTA file. Also returns a list
    of seqID's in the same order as that in the source FASTA file.
    """
    seqids = []
    seqid = None
    proteins = OrderedDict()
    for l in open(fasta):
        if l.startswith(">"):
            seqid, name = parse_fasta_header(l)
            seqids.append(seqid)
            proteins[seqid] = {
                'seq': "",
                'name': name,
                'original_header': l[1:].rstrip('\n').rstrip('\r'),
            }
            continue
        if seqid is not None:
            words = l.split()
            if words:
                proteins[seqid]['seq'] += words[0]
        proteins[seqid]['sequence_length'] = len(proteins[seqid]['seq'])

    proteins, id_mapping = generate_safe_seqids(proteins)

    return seqids, proteins


def print_proteins(proteins):
    """
    Prints out the proteins dictionary in a formatted
    manner that is also Python-eval compatible.
    """
    print "{"
    for seqid in proteins:
        print "  '%s': {" % seqid
        for key, value in proteins[seqid].items():
            print "    '%s': %s, " % (key, repr(value))
        print "  },"
    print "}"
    # Standard Library alternative
    # import pprint
    # pp = pprint.PrettyPrinter(indent=4)
    # print pp.pformat(proteins)


def write_proteins_fasta(
        fasta_filename, proteins, seqids, width=50):
    """
    Creates a fasta file of the sequences of a subset of the proteins.
    """
    f = open(fasta_filename, "w")
    out = proteins_to_fasta(proteins, seqids=seqids, width=width)
    f.write(out)
    f.close()


def proteins_to_fasta(proteins, seqids=[], use_safe_seqid=False, width=50):
    """
    Takes a proteins dictionary and returns a string containing
    all the sequences in FASTA format. Option parameters are
    a list of seqids to output (seqids) and the line width (width).
    """
    if seqids:
        idlist = seqids
    else:
        idlist = proteins

    fasta_out = ""
    for seqid in idlist:
        seq_wrap = textwrap.fill(proteins[seqid]['seq'], width)
        if use_safe_seqid:
            header = proteins[seqid]['safe_seqid']
        else:
            header = proteins[seqid]['name']
        fasta_out += ">%s\n%s\n" % (header, seq_wrap)
    return fasta_out


def chop_nterminal_peptide(protein, i_cut):
    """
    Finds all transmembrane fields in protein ("*_loops" and "*_helices")
    and deletes an N-terminal section of the protein indicated by i_cut.

    Assuming the topology of inner and outer loops remain the same, this
    may delete a certain number of elements at the N-terminus. Allows a
    simple way of removing lipo-protein and secretion-protein signals from
    the evaluation of the outer_loops of the protein.
    """
    protein['sequence_length'] -= i_cut
    for prop in protein:
        if '_loops' in prop or '_helices' in prop:
            sses = protein[prop]
            for i in range(len(sses)):
                j, k = sses[i]
                sses[i] = (j - i_cut, k - i_cut)
    for prop in protein:
        if '_loops' in prop or '_helices' in prop:
            sses = protein[prop]
            for i in reversed(range(len(sses))):
                j, k = sses[i]
                # tests if this loop or TM-helix has been cut out
                if j <= 0 and k <= 0:
                    del sses[i]
                # otherewise, neg value means loop or TM-helix is at the new N-terminal
                elif j <= 0 and k > 0:
                    sses[i] = (1, k)
                    # if a TM-helix gets cut in half and becomes a new N-terminal,
                    # convert the remaining residues to a loop and remove the helix
                    if '_helices' in prop:
                        program = prop.split('_')[0]
                        for x in protein:
                            if x == '%s_loops' % program:
                                new_N_loop = protein[x][0]
                                new_N_loop[0] = 1
                        del sses[i]


def generate_safe_seqids(proteins):
    """
    Takes a 'proteins' dictionary, keyed by sequence id, and
    adds a 'safe' sequence id ('safe_seqid') that is less likely
    to be munged or break various external programs.

    Returns a tuple of the updated proteins dictionary and a dictionary
    mapping the 'safe' sequence ids to the original ids.
    """
    id_mapping = {}
    count = 0
    for seqid in proteins:
        safe_id = re.sub(r'[^\w]', "", seqid) + '_' + `count`
        id_mapping[safe_id] = seqid
        proteins[seqid]['safe_seqid'] = safe_id
        count += 1

    return (proteins, id_mapping)


def clean_directory(top, excluded_files):
    """
    Deletes an entire directory tree, equivalent to rm -rf <top>,
    with an excluded file list.
    """
    for root, dirs, files in os.walk(top, topdown=False):
        for name in files:
            if name not in excluded_files:
                os.remove(os.path.join(root, name))
        for name in dirs:
            os.rmdir(os.path.join(root, name))


# Adapted from: https://stackoverflow.com/a/24618186/77990
def html2text(page, aggressive=False):
    soup = BeautifulSoup(page, 'html.parser')

    # kill all script and style elements
    for script in soup(["script", "style"]):
        script.extract()

    # get text
    text = soup.get_text()

    # break into lines and remove leading and trailing space on each
    lines = (line.strip() for line in text.splitlines())
    if aggressive:
        # break multi-headlines into a line each
        chunks = (phrase.strip() for line in lines for phrase in
                  line.split("  "))
        # drop blank lines
        text = '\n'.join(chunk for chunk in chunks if chunk)
    else:
        text = '\n'.join(lines)

    return text
