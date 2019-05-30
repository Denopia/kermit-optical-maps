import argparse
import os
from file_opener import open_valouev

'''
Colors contigs according to the output
of optical map mappin tool by Valouev et al.

This step could be merged with read coloring
so no "unnecessary" contig coloring file is 
generated.
'''


# Function to parse arguments given to this program
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='optical_mapping_output', action='store',
                        help='Path to the optical mapping output')
    parser.add_argument(dest='ref_opt_map_path', action='store',
                        help='Path to the reference optical maps')
    parser.add_argument(dest='contig_opt_map_path', action='store',
                        help='Path to the contig optical maps')
    parser.add_argument('-c', '--contig_colors_path', dest='contig_coloring_path', action='store', default="use_default",
                        help='Path where to save contig colors')
    parser.add_argument('-t', '--t_limit', action='store', dest='t_limit', default='0.0', type=float,
                        help='Minimum t-score')
    parser.add_argument('-s', '--s_limit', action='store', dest='s_limit', default='0.0', type=float,
                        help='Minimum s-score')
    results = parser.parse_args()
    return results


# Function to extract contigs and reference chromosome names from an output line
def extract_names(op_line):
    names = op_line.split("ment for ")[-1]
    contig_name = names.split(" and")[0].strip()
    reference_name = names.split("and ")[-1].strip()
    return contig_name, reference_name


# Function to extract contig fragment colors from an output line.
# Also returns updated values for the first and last mapping fragments.
def extract_colors(output_line, firstie, lastie):
    contig_side = output_line.split("]->")[0][1:].strip()
    reference_side = output_line.split("->[")[-1].replace("]", "").strip()
    contig_side_frags = contig_side.split(", ")
    reference_side_frags = reference_side.split(", ")

    temp_new_frags = []
    temp_new_frag_colors = []
    temp_firstie = firstie
    temp_lastie = lastie

    contig_side_total_length = 0.0
    for contig_frag in contig_side_frags:
        contig_side_total_length += float(contig_frag.split(":")[-1])
        if int(contig_frag.split(":")[0]) < temp_firstie:
            temp_firstie = int(contig_frag.split(":")[0])
        if int(contig_frag.split(":")[0]) > temp_lastie:
            temp_lastie = int(contig_frag.split(":")[0])

    reference_side_total_length = 0.0
    for reference_frag in reference_side_frags:
        reference_side_total_length += float(reference_frag.split(":")[-1])
        temp_new_frag_colors.append(int(reference_frag.split(":")[0]))
        temp_new_frags.append(float(reference_frag.split(":")[-1]))

    adjusted_new_frags = [(x*contig_side_total_length)/reference_side_total_length for x in temp_new_frags]

    return adjusted_new_frags, temp_new_frag_colors, temp_firstie, temp_lastie


# Function to finalize coloring of a single contig after its mapping information has been completely read
def finalize_coloring(contig_name, reference_name, new_contig, new_colors, full_contigs,
                      forward, firstie, lastie, reference_chromosome_nums, chromosome_shift):
    if forward:
        missing_from_start = firstie
        missing_from_end = len(full_contigs[contig_name]) - lastie - 1
        finalized_contigs = full_contigs[contig_name][:missing_from_start] + [x for x in new_contig]
        finalized_contigs = finalized_contigs + full_contigs[contig_name][len(full_contigs[contig_name])-missing_from_end:]
        #finalized_colors = missing_from_start * [new_colors[0] - 1] + [x for x in new_colors]
        finalized_colors = missing_from_start * [new_colors[0]] + [x for x in new_colors]
        #finalized_colors = finalized_colors + missing_from_end * [new_colors[-1] + 1]
        finalized_colors = finalized_colors + missing_from_end * [new_colors[-1]]

    if not forward:
        missing_from_end = firstie
        missing_from_start = len(full_contigs[contig_name]) - lastie - 1
        finalized_contigs = [x for x in new_contig]
        finalized_contigs.reverse()
        finalized_colors = [x for x in new_colors]
        finalized_colors.reverse()
        finalized_contigs = full_contigs[contig_name][:missing_from_start] + finalized_contigs
        finalized_contigs = finalized_contigs + full_contigs[contig_name][len(full_contigs[contig_name])-missing_from_end:]
        #finalized_colors = missing_from_start * [finalized_colors[0] + 1] + finalized_colors
        finalized_colors = missing_from_start * [finalized_colors[0]] + finalized_colors
        #finalized_colors = finalized_colors + missing_from_end * [finalized_colors[-1] - 1]
        finalized_colors = finalized_colors + missing_from_end * [finalized_colors[-1] - 1]

    reference_shift = reference_chromosome_nums[reference_name]
    shifted_finalized_colors = [reference_shift * chromosome_shift + int(x) for x in finalized_colors]

    return finalized_contigs, shifted_finalized_colors


def main(args):
    contig_reference_mapping = args.optical_mapping_output
    reference_opt_map_full = args.ref_opt_map_path
    contig_opt_map_full = args.contig_opt_map_path
    t_score_limit = args.t_limit
    s_score_limit = args.s_limit
    colored_contigs_path = args.contig_colors_path

    print("== Reading input files ==")
    # Read optical maps
    full_contigs = open_valouev(contig_opt_map_full, as_dict=True)
    full_reference = open_valouev(reference_opt_map_full, as_dict=True)

    # Read optical mapping output
    with open(contig_reference_mapping, "r") as rfile:
        mapping_lines = rfile.readlines()
    print("-- Done --")

    print("== Initializing ==")
    if colored_contigs_path == "use_default":
        if not os.path.exists("./data"):
            os.makedirs("./data")
        colored_contigs_path = "./data/colored_contigs.cval"

    # Find the right amount of padding between chromosomes

    # Dictionary: chromosome -> its number
    reference_chromosome_nums = {}
    # Start numbering chromosomes from 1
    rci = 1
    # Maximum number of fragments in any chromosome optical map
    max_fragments = 0
    for key in full_reference.keys():
        reference_chromosome_nums[key] = rci
        if len(full_reference[key]) > max_fragments:
            max_fragments = len(full_reference[key])
        rci += 1

    # Padding is the second smallest power of ten that is greater than the maximum number of chromosome fragments
    chromosome_shift = int("1" + (len(str(max_fragments)) + 1) * "0")

    # Store contig coloring information in these
    new_contig_names = []
    new_contig_fragments = []
    new_contig_fragments_colors = []

    # Mapping is forward/reverse
    forward = True
    # The first and last fragments of the mapping
    firstie = max_fragments
    lastie = 0
    # Coloring information of the current contig
    new_contig = []
    new_colors = []
    contig_name = "Name:NULL"
    reference_name = "Name:NULL"
    s_score = 0
    t_score = 0
    print("-- Done --")

    # Mapping output parsing loop
    print("== Parsing mapping output ==")
    for line in mapping_lines:
        # Next mapping is reversed
        if line.startswith("reverse"):
            forward = False
        # Next mapping is forward
        if line.startswith("forward"):
            forward = True
        # Line contains next mapping's contig and chromosome
        if line.startswith("alignment for"):
            contig_name, reference_name = extract_names(line)
        # Line contains fragment mapping information
        if line.startswith("["):
            additional_fragments, additional_fragment_colors, firstie, lastie = extract_colors(line, firstie, lastie)
            new_contig += additional_fragments
            new_colors += additional_fragment_colors
        # Line contains mapping's s-score
        if line.startswith("s-score"):
            s_score = float(line.split(":")[-1])
        # Line contains mapping's t-score, this line also ends the current mapping
        if line.startswith("t-score"):
            t_score = float(line.split(":")[-1])
            # If the mapping scores are good enough contig coloring is finalized and saved
            if t_score >= t_score_limit and s_score >= s_score_limit:
                final_fragments, final_colors = finalize_coloring(contig_name, reference_name, new_contig, new_colors,
                                                                  full_contigs, forward, firstie, lastie,
                                                                  reference_chromosome_nums, chromosome_shift)
                new_contig_names.append(contig_name)
                new_contig_fragments.append(final_fragments)
                new_contig_fragments_colors.append(final_colors)
            # Reset variables
            forward = True
            firstie = max_fragments
            lastie = 0
            new_contig = []
            new_colors = []
            contig_name = "Name:NULL"
            reference_name = "Name:NULL"
            s_score = 0
            t_score = 0
    print("-- Done --")

    # Write colored contigs in a file
    print("== Writing contig colors in a file ==")
    with open(colored_contigs_path, "w") as wfile:
        for i in range(len(new_contig_names)):
            name = new_contig_names[i]
            fragments = new_contig_fragments[i]
            colors = new_contig_fragments_colors[i]
            wfile.write(name + "\n")
            for frag in fragments:
                wfile.write(str(round(frag, 3)) + "\t")
            wfile.write("\n")
            for color in colors:
                wfile.write(str(color) + "\t")
            wfile.write("\n\n")
    print("-- Done --")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)

