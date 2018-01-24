import os
from inmembrane.helpers import dict_get, chop_nterminal_peptide


def get_annotations(params):
    """
    Creates a list of annotation functions required
    by this gram_pos protocol. The main program
    will run the annotation functions of this list,
    mapping the correct functions to the strings.

    As well, the function does some bookeeping on
    params to make sure the 'hmm_profiles_dir' is
    pointing in the right place.
    """
    annotations = []

    params['signalp4_organism'] = 'gram+'

    if not params['signalp4_bin'] or params[
        'signalp4_bin'] == 'signalp_scrape_web':
        annotations += ['signalp_scrape_web']
    else:
        annotations += ['signalp4']

    if not params['lipop1_bin'] or params['lipop1_bin'] == 'lipop_scrape_web':
        annotations += ['lipop_scrape_web']
    else:
        annotations += ['lipop1']

    annotations += ['hmmsearch3']

    if dict_get(params, 'helix_programs'):
        if 'tmhmm' in params['helix_programs']:
            if not params['tmhmm_bin'] or params[
                'tmhmm_bin'] == 'tmhmm_scrape_web':
                annotations.append('tmhmm_scrape_web')
            else:
                annotations.append('tmhmm')
        if 'memsat3' in params['helix_programs']:
            annotations.append('memsat3')

    params['hmm_profiles_dir'] = os.path.join(
        os.path.dirname(__file__), 'gram_pos_profiles')

    return annotations


def eval_surface_exposed_loop(
        sequence_length, n_transmembrane_region, outer_loops,
        terminal_exposed_loop_min, internal_exposed_loop_min):
    """
    This is the key algorithm in SurfG+ to identify Potentially Surface
    Exposed proteins. It evaluates all loops that poke out of the periplasmic
    side of the Gram+ bacterial membrane and tests if the loops are long
    enough to stick out of the peptidoglycan layer to be cleaved by proteases
    in a cell-shaving experiment.

    Returns True if any outer loop, or the N- or C-terminii are
    longer than the given thresholds.
    """

    outer_loops = outer_loops[:]

    if n_transmembrane_region == 0:
        # treat protein as one entire exposed loop
        return sequence_length >= terminal_exposed_loop_min

    if not outer_loops:
        return False

    loop_len = lambda loop: abs(loop[1] - loop[0]) + 1

    # if the N-terminal loop sticks outside
    if outer_loops[0][0] == 1:
        nterminal_loop = outer_loops[0]
        del outer_loops[0]
        if loop_len(nterminal_loop) >= terminal_exposed_loop_min:
            return True

    # if the C-terminal loop sticks outside
    if outer_loops:
        if outer_loops[-1][-1] == sequence_length:
            cterminal_loop = outer_loops[-1]
            del outer_loops[-1]
            if loop_len(cterminal_loop) >= terminal_exposed_loop_min:
                return True

    # test remaining outer loops for length
    for loop in outer_loops:
        if loop_len(loop) >= internal_exposed_loop_min:
            return True

    return False


def max_exposed_loop(
        sequence_length, n_transmembrane_region, outer_loops,
        terminal_exposed_loop_min, internal_exposed_loop_min):
    """
    This is the key algorithm in SurfG+ to identify Potentially Surface
    Exposed proteins. It evaluates all loops that poke out of the periplasmic
    side of the Gram+ bacterial membrane and tests if the loops are long
    enough to stick out of the peptidoglycan layer to be cleaved by proteases
    in a cell-shaving experiment.

    Returns True if any outer loop, or the N- or C-terminii are
    longer than the given thresholds.
    """

    outer_loops = outer_loops[:]

    if n_transmembrane_region == 0:
        # treat protein as one entire exposed loop
        return sequence_length

    if not outer_loops:
        return 0

    loop_len = lambda loop: abs(loop[1] - loop[0]) + 1

    lengths = []
    # if the N-terminal loop sticks outside
    if outer_loops[0][0] == 1:
        nterminal_loop = outer_loops[0]
        del outer_loops[0]
        lengths.append(loop_len(nterminal_loop))

    # if the C-terminal loop sticks outside
    if outer_loops:
        if outer_loops[-1][-1] == sequence_length:
            cterminal_loop = outer_loops[-1]
            del outer_loops[-1]
            lengths.append(loop_len(cterminal_loop))

    # test remaining outer loops for length
    for loop in outer_loops:
        lengths.append(loop_len(loop) / 2)

    return max(lengths)


def post_process_protein(params, protein):
    """
    This is the main analysis of the protein, where theprotein
    dictionary should contain all the necessary information
    from the annotations. Thus post_process_protein contain
    can determine the final analysis.
    """

    def sequence_length(protein):
        return protein['sequence_length']

    def has_tm_helix(protein):
        for program in params['helix_programs']:
            if dict_get(protein, '%s_helices' % program):
                return True
        return False

    def has_surface_exposed_loop(protein):
        for program in params['helix_programs']:
            if eval_surface_exposed_loop(
                    protein['sequence_length'],
                    len(protein['%s_helices' % (program)]),
                    protein['%s_outer_loops' % (program)],
                    params['terminal_exposed_loop_min'],
                    params['internal_exposed_loop_min']):
                return True
        return False

    def exposed_loop_extent(protein):
        extents = []
        for program in params['helix_programs']:
            if program + '_helices' in protein:
                extents.append(max_exposed_loop(
                    protein['sequence_length'],
                    len(protein['%s_helices' % (program)]),
                    protein['%s_outer_loops' % (program)],
                    params['terminal_exposed_loop_min'],
                    params['internal_exposed_loop_min']))
        if extents:
            return max(extents)
        else:
            return 0

    terminal_exposed_loop_min = \
        params['terminal_exposed_loop_min']

    is_hmm_profile_match = dict_get(protein, 'hmmsearch')
    is_lipop = dict_get(protein, 'is_lipop')
    if is_lipop:
        i_lipop_cut = protein['lipop_cleave_position']
    is_signalp = dict_get(protein, 'is_signalp')
    if is_signalp:
        i_signalp_cut = protein['signalp_cleave_position']

    details = []
    if is_hmm_profile_match:
        details += ["hmm(%s)" % "|".join(protein['hmmsearch'])]
    if is_lipop:
        details += ["lipop"]
    if is_signalp:
        details += ["signalp"]
    for program in params['helix_programs']:
        if has_tm_helix(protein):
            n = len(protein['%s_helices' % program])
            details += [program + "(%d)" % n]

    if is_lipop:
        chop_nterminal_peptide(protein, i_lipop_cut)
    elif is_signalp:
        chop_nterminal_peptide(protein, i_signalp_cut)

    if is_hmm_profile_match:
        category = "PSE-Cellwall"
    elif has_tm_helix(protein):
        if has_surface_exposed_loop(protein):
            category = "PSE-Membrane"
        else:
            category = "MEMBRANE(non-PSE)"
    else:
        if is_lipop:
            # whole protein considered outer terminal loop
            if sequence_length(protein) < terminal_exposed_loop_min:
                category = "LIPOPROTEIN(non-PSE)"
            else:
                category = "PSE-Lipoprotein"
        elif is_signalp:
            category = "SECRETED"
        else:
            category = "CYTOPLASM(non-PSE)"

    if details == []:
        details = ["."]

    protein['details'] = details
    protein['category'] = category
    if 'CYTOPLASM' not in category and 'SECRETED' not in category:
        protein['loop_extent'] = exposed_loop_extent(protein)
    else:
        protein['loop_extent'] = "."

    return details, category


def protein_output_line(seqid, proteins):
    return '%-15s   %-18s  %-4s %-50s  %s' % \
           (seqid,
            proteins[seqid]['category'],
            proteins[seqid]['loop_extent'],
            ";".join(proteins[seqid]['details']),
            proteins[seqid]['name'][:60])


def protein_csv_line(seqid, proteins):
    return '%s,%s,%s,%s,"%s"\n' % \
           (seqid,
            proteins[seqid]['category'],
            proteins[seqid]['loop_extent'],
            ";".join(proteins[seqid]['details']),
            proteins[seqid]['name'])


def summary_table(params, proteins):
    """
    Returns a string representing a simple summary table of
    protein classifcations.
    """
    out = ""
    counts = {}
    for seqid in proteins:
        category = proteins[seqid]['category']

        if category not in counts:
            counts[category] = 1
        else:
            counts[category] += 1

    out += "\n\n# Number of proteins in each class:\n"
    pse_total = 0

    for c in counts:
        if "PSE-" in c:
            pse_total += counts[c]
    counts["PSE(total)"] = pse_total

    for c in sorted(counts):
        if "PSE-" in c:
            spacer = "  "
            pse_total += counts[c]
        else:
            spacer = ""
        out += "# %s%-15s\t%i\n" % (spacer, c, counts[c])

    # out += "#\n# PSE(total)\t\t%i\n" % (pse_total)

    return out
