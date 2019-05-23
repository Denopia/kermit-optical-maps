import sys

read_contig_mapping = sys.argv[1]
contig_colors = sys.argv[2]
read_colors = sys.argv[3]
extend_colorings = sys.argv[4]
strict_coloring = sys.argv[5] #only first
alignment_limit = sys.argv[6] #min alignment length limit used

align_multiplier = 2

with open(contig_colors, "r") as rfile:
    contig_color_lines = rfile.readlines()

with open(read_contig_mapping, "r") as rfile:
    mapping_lines = rfile.readlines()

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

#print(contig_fragments_colors)

read_to_contig = {}
qq = 0
last_read_name = "Read:NULL"
for line in mapping_lines:
    qq += 1
    mapping = line.split("\t")
    if qq % 1000 == 0:
        print(line)
    '''
    print("Read name:", mapping[0],  "Read length: ", mapping[1], "Read start:", mapping[2], "Read end:", mapping[3],
          "Contig name:", mapping[5], "Contig length:", mapping[6], "Contig start:", mapping[7], "Contig end:", mapping[8],
          "Matches:", mapping[9], "Match/mismatch/indel:", mapping[10])
    print()
    '''
    if contig_fragments_colors.get(mapping[5], -1) != -1 and \
            (strict_coloring == "0" or (strict_coloring == "1" and mapping[0] != last_read_name)):

        if read_to_contig.get(mapping[0], -1) == -1:
            if (alignment_limit == "1" and align_multiplier * int(mapping[10]) >= int(mapping[1])) or (alignment_limit == "0"):
                read_to_contig[mapping[0]] = mapping[:11]
        elif int(read_to_contig[mapping[0]][9]) < int(mapping[9]):
            if (alignment_limit == "1" and align_multiplier * int(mapping[10]) >= int(mapping[1])) or (alignment_limit == "0"):
                read_to_contig[mapping[0]] = mapping[:11]

    last_read_name = mapping[0]

read_keys = read_to_contig.keys()
new_keys = []
for key in read_keys:
    new_keys.append(int(key.split("_")[-1].strip()))

new_keys.sort()

newest_keys = []
for key in new_keys:
    newest_keys.append("Read_" + str(key))

reads_and_colors = []

for key in newest_keys:
    mapping = read_to_contig[key]
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

    if extend_colorings == "0":
        start = contig_start
        end = contig_end
    else:
        start = contig_start - read_start
        end = contig_end + (read_length - read_end)

    current_length = 0
    current_fragment = 0
    start_color = 999999999
    end_color = -1

    while current_fragment < len(contig_fragments):
        current_length += 1000 * contig_fragments[current_fragment]
        if current_length > start and start_color == 999999999:
            start_color = contig_colors[current_fragment]
        if current_length > end and end_color == -1:
            end_color = contig_colors[current_fragment]
        current_fragment += 1

    if end_color == -1:
        end_color = contig_colors[-1]

    if start_color > end_color:
        temp_start = start_color
        start_color = end_color
        end_color = temp_start


    reads_and_colors.append([read_name, start_color, end_color])


with open(read_colors, "w") as wfile:
    for coloring in reads_and_colors:
        wfile.write(coloring[0] + "\t" + str(coloring[1]) + "\t" + str(coloring[2]) + "\n")
