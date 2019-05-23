import sys
import argparse
from file_opener_f import open_valouev

# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument(dest='optical_mapping_result', action='store',
                    help='Path to optical mapping output')

parser.add_argument(dest='optmap_path', action='store',
                    help='Path where to save optical maps')

parser.add_argument('-z', '--enzyme', action='store', dest='enzyme', default='XhoI',
                    help='Simulated restriction enzyme (str)')

parser.add_argument('-f', '--min_fragment', action='store', dest='min_fragment',
                    help='Minimum fragment length (int)', default='50', type=int)

parser.add_argument('-b', '--skip_first', action='store_true', dest='skip_first',
                    help='Ignore first fragment', default=False)

parser.add_argument('-e', '--skip_last', action='store_true', dest='skip_last',
                    help='Ignore last fragment', default=False)

parser.add_argument('-t', '--t_limit', action='store', dest='t_limit',
                    help='Minimum t-score', default="0.0", type=float)

parser.add_argument('-s', '--s_limit', action='store', dest='s_limit',
                    help='Minimum s-score', default="0.0", type=float)



results = parser.parse_args()


# READ ARGUMENTS
contig_reference_mapping = sys.argv[1]

contig_opt_map_trim = sys.argv[2]
reference_opt_map_trim = sys.argv[3]
contig_opt_map_full = sys.argv[4]
reference_opt_map_full = sys.argv[5]
colored_contigs = sys.argv[6]

t_score_limit = int(sys.argv[7])
s_score_limit = int(sys.argv[8])

# OPEN FILES
trimmed_contigs = open_valouev(contig_opt_map_trim, as_dict=True)
trimmed_reference = open_valouev(reference_opt_map_trim, as_dict=True)
full_contigs = open_valouev(contig_opt_map_full, as_dict=True)
full_reference = open_valouev(reference_opt_map_full, as_dict=True)

# SAVE COLORED FRAGMENT INFORMATION IN THESE
new_contig_names = []
new_contig_fragments = []
new_contig_fragments_colors = []

#t_score_limit = 15
#s_score_limit = 30

#t_score_limit = 0
#s_score_limit = 0


# OPEN MAPPING OUTPUT
with open(contig_reference_mapping, "r") as rfile:
    mapping_lines = rfile.readlines()

# FIND THE CORRECT AMOUNT OF DISTANCE TO ADD BETWEEN COLORINGS IN DIFFERENT CHROMOSOMES
reference_chromosome_nums = {}
rci = 1
max_fragments = 0
for key in full_reference.keys():
    reference_chromosome_nums[key] = rci
    if len(full_reference[key]) > max_fragments:
        max_fragments = len(full_reference[key])
    rci += 1

chromosome_shift = int("1" + (len(str(max_fragments)) + 1) * "0")

# FUNCTION THAT EXTRACTS CONTIG + REFERENCE CHROMOSOME NAMES FROM A LINE
def extract_names(line):
    names = line.split("ment for ")[-1]
    contig_name = names.split(" and")[0].strip()
    reference_name = names.split("and ")[-1].strip()
    return contig_name, reference_name


# FUNCTION THAT EXTRACTS CONTIG FRAGMENT COLORS FROM A LINE
def extract_colors(line, frags, frag_colors, firstie, lastie):
    contig_side = line.split("]->")[0][1:].strip()
    reference_side = line.split("->[")[-1].replace("]", "").strip()
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

    return frags + adjusted_new_frags, frag_colors + temp_new_frag_colors, temp_firstie, temp_lastie


# FUNCTION THAT FINALIZES COLOR INFORMATION OF A SINGLE CONTIG AFTER WHOLE MAPPING IS READ
def finalize_coloring(contig_name, reference_name, new_contig, new_colors,
                          trimmed_contigs, full_contigs, forward, firstie, lastie):

    if forward == -1:
        print("!!!!!!!!!!!!!!!!!!")
        print("!!!            !!!")
        print("!!! HUGE ERROR !!!")
        print("!!!            !!!")
        print("!!!!!!!!!!!!!!!!!!")
        quit()

    if forward == 1:
        missing_from_start = 1 + firstie
        missing_from_end = len(full_contigs[contig_name]) - 1 - (lastie + 1)

        finalized_contigs = full_contigs[contig_name][:missing_from_start] + [x for x in new_contig]
        finalized_contigs = finalized_contigs + full_contigs[contig_name][len(full_contigs[contig_name])-missing_from_end:]

        finalized_colors = missing_from_start * [new_colors[0] - 1] + [x for x in new_colors]
        finalized_colors = finalized_colors + missing_from_end * [new_colors[-1] + 1]

    if forward == 0:
        missing_from_end = 1 + firstie
        missing_from_start = len(full_contigs[contig_name]) - 1 - (lastie + 1)
        finalized_contigs = [x for x in new_contig]
        finalized_contigs.reverse()
        finalized_colors = [x for x in new_colors]
        finalized_colors.reverse()
        finalized_contigs = full_contigs[contig_name][:missing_from_start] + finalized_contigs
        finalized_contigs = finalized_contigs + full_contigs[contig_name][len(full_contigs[contig_name])-missing_from_end::]

        finalized_colors = missing_from_start * [finalized_colors[0] + 1] + finalized_colors
        finalized_colors = finalized_colors + missing_from_end * [finalized_colors[-1] - 1]

    reference_shift = reference_chromosome_nums[reference_name]
    shifted_finalized_colors = [reference_shift * chromosome_shift + int(x) for x in finalized_colors]

    return finalized_contigs, shifted_finalized_colors


# SOME VARIABLES NEEDED FOR MAPPING PARSING
forward = -1
firstie = 999999
lastie = 0
new_contig = []
new_colors = []
contig_name = "Name:NULL"
reference_name = "Name:NULL"


# MAPPING PARSING LOOP
for line in mapping_lines:
    if line.startswith("reverse"):
        forward = 0
    if line.startswith("forward"):
        forward = 1
    if line.startswith("alignment for"):
        contig_name, reference_name = extract_names(line)
    if line.startswith("["):
        new_contig, new_colors, firstie, lastie = extract_colors(line, new_contig, new_colors, firstie, lastie)
    if line.startswith("s-score"):
        s_score = float(line.split(":")[-1])
    if line.startswith("t-score"):
        t_score = float(line.split(":")[-1])
        if t_score > t_score_limit and s_score > s_score_limit:
            final_fragments, final_colors= finalize_coloring(contig_name, reference_name, new_contig, new_colors,
                                                      trimmed_contigs, full_contigs, forward, firstie, lastie)
            new_contig_names.append(contig_name)
            new_contig_fragments.append(final_fragments)
            new_contig_fragments_colors.append(final_colors)
        forward = -1
        firstie = 999999
        lastie = 0
        new_contig = []
        new_colors = []
        contig_name = "Name:NULL"
        reference_name = "Name:NULL"
        s_score = 0
        t_score = 0


# WRITE FINALIZED COLORING INFORMATION IN A FILE
with open(colored_contigs, "w") as wfile:
    for i in range(len(new_contig_names)):
        name = new_contig_names[i]
        fragments = new_contig_fragments[i]
        colors = new_contig_fragments_colors[i]
        wfile.write(name + "\n")
        for frag in fragments:
            wfile.write(str(round(frag, 3)) + "\t")
        wfile.write("\n")
        # ADD 1 TO SHIFT COLORS BECAUSE THEY START FROM 0 IN THE TRIMMED VERSION OF REFERENCE
        for color in colors:
            wfile.write(str(color+1) + "\t")
        wfile.write("\n\n")







