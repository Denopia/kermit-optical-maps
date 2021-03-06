'''

Copyright 2018, 2019 Miika Leinonen

This file is part of OpticalKermit.

OpticalKermit is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OpticalKermit is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with OpticalKermit.  If not, see <https://www.gnu.org/licenses/>.

'''

import argparse
import os
from file_opener import open_cf


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='read_colorings', action='store',
                        help='Path to the read colors')
    parser.add_argument('-a', '--adjusted_colors_path', dest='adj_colors', action='store', default="use_default",
                        help='Path where to save adjusted read colors')
    results = parser.parse_args()
    return results


def find_index(color, clist):
    for i in range(len(clist)):
        if clist[i] == color:
            return i
    return 0


def main(args):

    print("=== Read recoloring started ===")
    print("== Initializing ==")
    coloring_path = args.read_colorings
    new_coloring_path = args.adj_colors
    
    if new_coloring_path == "use_default":
        if not os.path.exists("./data"):
            os.makedirs("./data")
        new_coloring_path = "./data/colored_reads_adjusted.cf"

    colorings = open_cf(coloring_path)
    multiplier = int("1" + (len(str(colorings[0][1])) - 1) * "0")
    big_shift = len(str(multiplier)) - 1

    colors = {}
    #oo = 0
    print("-- Done --\n")

    print("== Finding used colors ==")
    for cr in colorings:
        #oo += 1
        chromos = str(cr[1])[:-big_shift]
        #print(oo, cr, chromos)
        if colors.get(chromos, "na") == "na":
            colors[chromos] = []
        first_color = cr[1] - int(chromos) * multiplier
        second_color = cr[2] - int(chromos) * multiplier
        if first_color not in colors[chromos]:
            colors[chromos].append(first_color)
        if second_color not in colors[chromos]:
            colors[chromos].append(second_color)

    for key in colors.keys():
        colors[key].sort()
    print("-- Done --\n")

    print("== Writing adjusted colorings ==")
    with open(new_coloring_path, 'w') as wfile:
        for cr in colorings:
            chromos = str(cr[1])[:-big_shift]
            first_color = cr[1] - int(chromos) * multiplier
            second_color = cr[2] - int(chromos) * multiplier

            si = find_index(first_color, colors[chromos])
            ei = find_index(second_color, colors[chromos])

            wfile.write(cr[0] + "\t" + str(si + int(chromos) * multiplier) + "\t" + str(ei + int(chromos) * multiplier) + "\n")

    print("-- Done --\n")
    print("=== Read recoloring finished ===\n")


if __name__ == "__main__":
    args = parse_arguments()
    main(args)

