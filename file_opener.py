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

#import pysam
import Bio.SeqIO
import Bio.Restriction


def open_fasta(path):
    fastafile = open(path, "rU")
    records = []
    for record in Bio.SeqIO.parse(fastafile, "fasta"):
        records.append(record)
    return records


def open_fastq(path):
    fastqfile = open(path, "rU")
    records = []
    for record in Bio.SeqIO.parse(fastqfile, "fastq"):
        records.append(record)
    return records


# Reads fragment file (.valouev format)
def open_valouev(path, include_name=False, as_dict=False):
    with open(path, 'r') as rfile:
        lines = rfile.readlines()

    if as_dict:
        fragment_dict = {}
        for line in lines:
            if len(line) < 2:
                continue
            if len(line.split("\t")) == 1:
                name = line[:-1]
                continue
            fragments = line.split("\t")
            frags = [float(x) for x in fragments[2:]]
            fragment_dict[name] = frags
        return fragment_dict

    else:
        fragment_lists = []
        for line in lines:
            if len(line) < 2:
                continue
            if len(line.split("\t")) == 1:
                name = line[:-1]
                continue
            fragments = line.split("\t")
            if include_name:
                frags = [name] + [float(x) for x in fragments[2:]]
            else:
                frags = [float(x) for x in fragments[2:]]
            fragment_lists.append(frags)
        return fragment_lists


# Reads read color file (.cf format)
def open_cf(path):
    with open(path, 'r') as rfile:
        colorings = rfile.readlines()
    coloring_list = []
    for line in colorings:
        rcp = line.split("\t")
        coloring_list.append([rcp[0], int(rcp[1]), int(rcp[2])])
    return coloring_list


"""
def open_sam(path):
    samfile = pysam.AlignmentFile(path, "rb")
    alignments = []
    for alignment in samfile.fetch():
        alignments.append(alignment)
    return alignments


def open_bam(path):
    bamfile = pysam.AlignmentFile(path, "rb", check_sq=False)
    alignments = []
    for alignment in bamfile.fetch(until_eof=True):
        alignments.append(alignment)
    return alignments
"""
