import argparse
import os

'''
Colors reads according to reads to contigs
alignments and contig colors.

Could probably be merged with contig coloring.
'''


# Function to parse arguments given to this program
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='read_contig_mapping', action='store',
                        help='Path to the reads to contigs alignments')
    parser.add_argument(dest='contig_colors', action='store',
                        help='Path to the contig colorings')
    parser.add_argument('-r', '--read_colors_path', dest='read_colors', action='store', default="use_default",
                        help='Path where to save read colors')
    parser.add_argument('-e', '--extend', action='store_true', dest='extend_coloring', default=False,
                        help='Extend coloring')
    parser.add_argument('-s', '--strict', action='store_true', dest='strict_coloring', default=False,
                        help='Color based on only the best alignment')
    parser.add_argument('-l', '--alignment_limit', action='store', dest='alignment_limit', default='0.0', type=float,
                        help='Alignment limit multiplier')
    results = parser.parse_args()
    return results


def main(args):
    print("=== Read coloring started ===")
    print("== Initializing ==")
    read_contig_mapping = args.read_contig_mapping
    contig_colors = args.contig_colors
    read_colors_path = args.read_colors
    extend_colorings = args.extend_coloring
    strict_coloring = args.strict_coloring
    alignment_limit_multiplier = args.alignment_limit

    if read_colors_path == "use_default":
        if not os.path.exists("./data"):
            os.makedirs("./data")
        read_colors_path = "./data/colored_reads.cf"

    print("-- Done --\n")

    print("== Reading input files ==")
    with open(contig_colors, "r") as rfile:
        contig_color_lines = rfile.readlines()

    with open(read_contig_mapping, "r") as rfile:
        mapping_lines = rfile.readlines()


    # Dictionary: Contig -> [[Fragment lengths], [Fragment colors]]
    contig_fragments_colors = {}
    i = 0
    while i < len(contig_color_lines):
        contig_name = contig_color_lines[i].strip()
        contig_fragments_colors[contig_name] = []
        i += 1
        fragments = contig_color_lines[i].strip().split("\t")
        i += 1
        colors = contig_color_lines[i].strip().split("\t")
        i += 2
        contig_fragments_colors[contig_name].append([float(x) for x in fragments])
        contig_fragments_colors[contig_name].append([int(x) for x in colors])
    print("-- Done --\n")

    print("== Coloring contigs ==")
    # Dictionary: Read -> Its chosen alignment to a contig
    read_to_contig = {}
    last_read_name = "Name:NULL"
    for line in mapping_lines:
        mapping = line.split("\t")
        '''
        print("Read name:", mapping[0],  "Read length: ", mapping[1], "Read start:", mapping[2],
              "Read end:", mapping[3], "Contig name:", mapping[5], "Contig length:", mapping[6],
              "Contig start:", mapping[7], "Contig end:", mapping[8], "Matches:", mapping[9],
              "Match/mismatch/indel:", mapping[10])
        '''
        if contig_fragments_colors.get(mapping[5], -1) != -1 and \
                (not strict_coloring or (strict_coloring and mapping[0] != last_read_name)):
            if read_to_contig.get(mapping[0], -1) == -1:
                if int(mapping[10]) >= alignment_limit_multiplier * int(mapping[1]):
                    read_to_contig[mapping[0]] = mapping[:11]
            elif int(read_to_contig[mapping[0]][10]) < int(mapping[10]):
                if int(mapping[10]) >= alignment_limit_multiplier * int(mapping[1]):
                    read_to_contig[mapping[0]] = mapping[:11]
        last_read_name = mapping[0]

    reads_and_colors = []

    # Color reads
    for read_name in sorted(read_to_contig.keys(), key=lambda item: (len(item), item)):
        mapping = read_to_contig[read_name]
        read_name = mapping[0]
        contig_name = mapping[5]
        read_length = int(mapping[1])
        read_start = int(mapping[2])
        read_end = int(mapping[3])
        contig_length = int(mapping[6])
        contig_start = int(mapping[7])
        contig_end = int(mapping[8])
        sign = mapping[4]

        contig_fragments = contig_fragments_colors[contig_name][0]
        contig_colors = contig_fragments_colors[contig_name][1]

        if not extend_colorings:
            start = contig_start
            end = contig_end
        else:
            start = contig_start - read_start
            end = contig_end + (read_length - read_end)

        current_length = 0
        current_fragment = 0
        start_color = -1
        end_color = -1

        while current_fragment < len(contig_fragments):
            current_length += 1000 * contig_fragments[current_fragment]
            if current_length > start and start_color == -1:
                start_color = contig_colors[current_fragment]
            if current_length > end and end_color == -1:
                end_color = contig_colors[current_fragment]
            current_fragment += 1

        # If the end color stays at -1 in some extreme cases set it to the last color of the contig
        if end_color == -1:
            end_color = contig_colors[-1]

        # If start color stays at -1 there is most likely some error and the coloring should not be saved
        if start_color == -1:
            continue

        # If the start and end color order is inverted flip it
        # (Due to reverse mapping in contig coloring)
        if start_color > end_color:
            temp_start = start_color
            start_color = end_color
            end_color = temp_start

        reads_and_colors.append([read_name, start_color, end_color])
    print("-- Done --\n")

    print("== Writing read colors in a file ==")
    with open(read_colors_path, "w") as wfile:
        for coloring in reads_and_colors:
            wfile.write(coloring[0] + "\t" + str(coloring[1]) + "\t" + str(coloring[2]) + "\n")
    print("-- Done --\n")
    print("=== Read coloring finished ===\n")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)

