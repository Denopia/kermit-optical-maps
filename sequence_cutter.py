import Bio.SeqIO
import Bio.Restriction


# Function to cut a single sequence into fragments.
# Supports single or multiple enzymes.
# Multiple enzymes in "Enzyme1_Enzyme2_Enzyme3_Enzyme4" format.
def fragmentize_sequence(sequence, enzymes="XhoI_NheI_EagI"):
    enzyme_batch = enzymes.split("_")
    restriction_batch = Bio.Restriction.RestrictionBatch(enzyme_batch)
    cut_sites = restriction_batch.search(sequence)
    merged_cut_sites = []
    for enz, cps in cut_sites.items():
        merged_cut_sites = merged_cut_sites + cps
    merged_cut_sites.sort()

    if len(merged_cut_sites) == 0:
        return []

    fragment_lengths = [merged_cut_sites[0]]
    for i in range(1, len(merged_cut_sites)):
        fragment_lengths.append(merged_cut_sites[i] - merged_cut_sites[i-1])
    fragment_lengths.append(len(sequence) - merged_cut_sites[-1])
    return fragment_lengths
