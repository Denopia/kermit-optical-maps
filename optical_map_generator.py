import argparse
from file_opener import open_fasta
from sequence_cutter import fragmentize_sequence


def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(dest='fasta_path', action='store',
                        help='Path to sequences')
    parser.add_argument('-o', '--optical_map_path', dest='optmap_path', action='store', default="use_default",
                        help='Path where to save optical maps')
    parser.add_argument('-z', '--enzyme', action='store', dest='enzyme', default='XhoI',
                        help='Simulated restriction enzyme (str)')
    parser.add_argument('-f', '--min_fragment', action='store', dest='min_fragment', default='50', type=int,
                        help='Minimum fragment length (int)')
    parser.add_argument('-b', '--skip_first', action='store_true', dest='skip_first', default=False,
                        help='Ignore first fragment')
    parser.add_argument('-e', '--skip_last', action='store_true', dest='skip_last', default=False,
                        help='Ignore last fragment')
    results = parser.parse_args()
    return results


def main(args):
    print("=== Optical map generation started ===")
    sequences_fasta_path = args.fasta_path
    optical_map_path = args.optmap_path
    enzyme = args.enzyme
    minimum_fragment_length = args.min_fragment
    skip_first = args.skip_first
    skip_last = args.skip_last
    if optical_map_path == "use_default":
        optical_map_path = sequences_fasta_path + ".valouev"

    # Open sequences
    print("== Reading sequences ==")
    sequences = open_fasta(sequences_fasta_path)
    print("-- Done --\n")

    # Fragmentize sequences
    print("== Fragmentizing sequences ==")
    sequences_fragments = {}
    for sequence in sequences:
        seq_name = sequence.id
        seq_fragments = fragmentize_sequence(sequence.seq, enzyme)
        if len(seq_fragments) > 0:
            sequences_fragments[seq_name] = seq_fragments
    print("-- Done --")

    # Write optical maps file
    print("== Writing optical maps ==")
    with open(optical_map_path, "w") as wfile:
        for sequence_name in sorted(sequences_fragments.keys(), key=lambda item: (len(item), item)):
            sequence_fragments = sequences_fragments[sequence_name]

            start_fragment = 0
            if skip_first:
                start_fragment += 1

            end_fragment = len(sequence_fragments) - 1
            if skip_last:
                end_fragment -= 1

            wfile.write(sequence_name + "\n" + enzyme + "\t" + enzyme)
            for i in range(start_fragment, end_fragment + 1):
                current_fragment = sequence_fragments[i]
                if current_fragment >= minimum_fragment_length:
                    wfile.write("\t" + str(current_fragment / 1000))

            wfile.write("\n\n")
    print("-- Done --\n")
    print("=== Optical map generation finished ===\n")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
