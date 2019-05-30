# Optical map extension to Kermit

This is an extension to the guided genome assembler Kermit https://github.com/rikuu/kermit. Kermit uses a reference genome or lonkage maps to color reads before assembling the genome. This extension uses optical maps to color reads, which can then be given to Kermit as input.

### Assembly steps
1. Reads are corrected (e.g. using correction tool CONSENT https://github.com/morispi/CONSENT)
2. The corrected reads are used to assemble pre-coloring contigs with miniasm https://github.com/lh3/miniasm
3. Optical maps of the pre-coloring contigs are generated using optical_map_generator.py (imitating restriction enzyme used for reference optical map, 'XhoI' as default)
4. Optical map for the reference is generated like in step 3 if it is not available but the genome itself is i.e. in test cases
5. Pre-coloring contig optical maps are mapped to the reference optical map using Valouev mapping tool https://github.com/mmuggli/valouev_optmap_alignment
6. Cotigs are colored using contig_colorer.py
7. Corrected reads are colored using read_colorer.py
8. Post-coloring contigs are assembled with Kermit using the colored reads

Alternatively you could try skipping correcting reads and polish the post-coloring reads at the end.

### External programs needed
1. Long read correction tool, for example CONSENT https://github.com/morispi/CONSENT
2. Minimap https://github.com/lh3/minimap2
3. Miniasm https://github.com/lh3/miniasm
4. Valouev mapping tool https://github.com/mmuggli/valouev_optmap_alignment
5. Kermit https://github.com/rikuu/kermit

### Other requirements
1. Python 3
2. Python libraries: Biopython https://biopython.org/
3. Anything needed for the external programs to work

### Installation
Through this guide we assume the home directory is ??? and everything is installed there
Where to install,
How to install (external apps)


### Usage
This extension contains three main programs, optical map generator (optical_map_generator.py), contig colorer (contig_colorer.py) and read colorer (read_colorer.py).

Optical map generator is used with 
```
python3 optical_map_generator.py [sequenses_path.fasta] [optional arguments] 
```
Required arguments:
```
sequences_path.fasta : Path to the sequences you want to generate optical map for.
```
Optional arguments:
```
-o optical_map_path.valouev : Path to where optical map file is saved, by default 'sequences_path.fasta.valouev'
-z enzyme : Name of the enzyme to be emulated, by default 'XhoI'. Multiple enzymes separated with '_' e.g. 'XhoI_NheI'
-f min_fragment_size : Minimum fragment size, by default 50
-b : Skips first fragment, included by default
-e : Skips last fragment, included by default
```

Read colorer is used with 
```
python ???.py a b c d
```
a explanation
b explanation
c explanation
d explanation

### Example run



