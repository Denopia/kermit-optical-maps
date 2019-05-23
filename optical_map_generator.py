import argparse
#from pathlib import Path
#import os
from file_opener import open_fasta
from sequence_cutter import fragmentize_sequence

# Parse arguments
parser = argparse.ArgumentParser()

parser.add_argument(dest='fasta_path', action='store',
                    help='Path to sequences')

parser.add_argument('-o', '--optical_map_path', dest='optmap_path', action='store', default="use_default",
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
if optical_map_path == "use_default":
    optical_map_path = sequences_fasta_path + ".valouev"


# Open sequences
print("== Reading sequences ==")
sequences = open_fasta(sequences_fasta_path)
print("-- Done --")

# Fragmentize sequences
print("== Fragmentizing sequences ==")
sequences_fragments = {}
for sequence in sequences:
    seq_name = sequence.id
    seq_fragments = fragmentize_sequence(sequence.seq, enzyme)
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
print("-- Done --")
