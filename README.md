# Optical map extension to Kermit

This is an extension to the guided genome assembler Kermit https://github.com/rikuu/kermit. Kermit uses a reference genome or lonkage maps to color reads before assembling the genome. This extension uses optical maps to color reads, which can then be given to Kermit as input.

## Assembly steps
1. Reads are corrected (e.g. using correction tool CONSENT https://github.com/morispi/CONSENT)
2. The corrected reads are used to assemble pre-coloring contigs with miniasm https://github.com/lh3/miniasm
3. Optical maps of the pre-coloring contigs are generated using optical_map_generator.py (imitating restriction enzyme used for reference optical map, 'XhoI' as default)
4. Optical map for the reference is generated like in step 3 if it is not available but the genome itself is i.e. in test cases
5. Pre-coloring contig optical maps are mapped to the reference optical map using Valouev mapping tool https://github.com/mmuggli/valouev_optmap_alignment
6. Cotigs are colored using contig_colorer.py
7. Corrected reads are colored using read_colorer.py
8. Read colors are adjusted to avoid gaps due to unused colors (optional)
9. Post-coloring contigs are assembled with Kermit using the colored reads

Alternatively you could try skipping correcting reads and polish the post-coloring reads at the end.

## External programs needed
1. Long read correction tool, for example CONSENT https://github.com/morispi/CONSENT
2. Minimap https://github.com/lh3/minimap2
3. Miniasm https://github.com/lh3/miniasm
4. Valouev mapping tool https://github.com/mmuggli/valouev_optmap_alignment
5. Kermit https://github.com/rikuu/kermit

### Other requirements
1. Python 3
2. Python libraries: Biopython https://biopython.org/
3. Anything needed for the external programs to work

## Installation
Installation guide might appear later. This extension has only executable python files, so you basically just need to install the external programs, python3 and Biopython according to their respective guides.


## Usage
This extension contains four main programs, optical map generator (optical_map_generator.py), contig colorer (contig_colorer.py), read colorer (read_colorer.py) and read recolorer (read_recolorer.py).

### Optical map generator 
```
python3 optical_map_generator.py [sequenses_path.fasta] [optional arguments] 
```
Required arguments:
```
sequences_path.fasta : Path to the sequences you want to generate optical map for, in fasta format.
```
Optional arguments:
```
-o optical_map_path.valouev : Path to where optical map file is saved to, by default 'sequences_path.fasta.valouev'
-z enzyme : Name of the enzyme to be emulated, by default 'XhoI'. Multiple enzymes separated with '_' e.g. 'XhoI_NheI'
-f min_fragment_size : Minimum fragment size, by default 50
-b : Skips first fragment, included by default
-e : Skips last fragment, included by default
```

### Contig colorer
```
python3 contig_colorer.py [optical_mapping_output.txt, reference_map_path.valouev, contigs_maps_path.valouev]
                          [optional arguments] 
```
Required arguments:
```
optical_mapping_output.txt : Path to the output of optical map mapping tool by Valouev et al.
reference_map_path.valouev : Path to the reference genome optical map, in .valouev format
contigs_maps_path.valouev : Path to the pre-coloring contigs optical maps, in .valouev format
```
Optional arguments:
```
-c colored_contigs_path.cval : Path to where colored contigs file is saved to, by default './data/colored_contigs.cval'
-t t_limit : t-score threshold that must be exceeded for the contig to be colored, by default 0
-s s_limit : s-score threshold that must be exceeded for the contig to be colored, by default 0
```

### Read colorer
```
python3 read_colorer.py [alignments_path.txt, colored_contigs_path.cval]
                        [optional arguments] 
```
Required arguments:
```
alignments_path.txt : Path to the reads-to-contigs alignment output of minimap2
colored_contigs_path.cval : Path to the colored contigs file, in .cval format
```
Optional arguments:
```
-r colored_reads_path.cf : Path to where colored reads file is saved to, by default './data/colored_reads.cf'
-e : extends read coloring, does not do this by default
-s : uses strict coloring i.e. coloring is done only based on the best aligning contig
-l len_multiplier : length of the aligning section must be at least len_multiplier * read_length for the read to be colored, by default 0.0
```


### Read recolorer
```
python3 read_recolorer.py [colored_reads_path.cf] [optional arguments] 
```
Required arguments:
```
colored_reads_path.cf : Path to the colored reads, in .cf format
```
Optional arguments:
```
-a adjusted_colored_reads_path.cf : Path to where adjusted colored reads file is saved to, by default './data/colored_reads_adjusted.cf'
```

### Example run
genome-assembly.sh contains an example run of the whole pipeline. Example run with explanations might appear in this readme later.

