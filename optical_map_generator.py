import argparse
from file_opener_f import open_fasta
from sequence_cutter_f import fragmentize_sequence

# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument(dest='fasta_path', action='store',
                    help='Path to sequences')

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

results = parser.parse_args()

#print(results)

sequences_fasta_path = results.fasta_path
optical_map_path = results.optmap_path
enzyme = results.enzyme
minimum_fragment_length = results.min_fragment
skip_first = results.skip_first
skip_last = results.skip_last


# Open sequences
sequences = open_fasta(sequences_fasta_path)

# Write fragments in a file
with open(optical_map_path, "w") as wfile:
    for sequence in sequences:
        sequence_name = sequence.id
        sequence_fragments = fragmentize_sequence(sequence.seq, enzyme)

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
