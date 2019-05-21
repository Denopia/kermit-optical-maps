# Optical map extension to Kermit

This is an extension to the guided genome assembler Kermit https://github.com/rikuu/kermit. Kermit uses a reference genome or lonkage maps to color reads before assembling the genome. This extension uses optical maps to color reads, which can then be given to Kermit as input.

### Here are the steps needed to color reads with optical maps:
1. Correct reads (we used correction tool CONSENT https://github.com/morispi/CONSENT)
2. The corrected reads are used to assemble pre-coloring contigs with miniasm https://github.com/lh3/miniasm
3.1. Optical maps of the pre-coloring contigs are generated using ???.py (imitating restriction enzyme used for reference optical map, 'XhoI' as default)
3.2. Optical map for the reference is generated like in step 3.1. if it is not available but the genome itself is i.e. in test cases
4. Pre-coloring contig optical maps are mapped to the reference optical map using Valouev mapping tool https://github.com/mmuggli/valouev_optmap_alignment
5. Corrected reads are colored using ???.py
6. Post-coloring contigs are assembled with Kermit using the colored reads


### External programs needed:
1. Long read correction tool, for example CONSENT https://github.com/morispi/CONSENT
2. Minimap https://github.com/lh3/minimap2
3. Miniasm https://github.com/lh3/miniasm
4. Valouev mapping tool https://github.com/mmuggli/valouev_optmap_alignment
5. Kermit https://github.com/rikuu/kermit
(Plus other requirements for these programs to work)
6. Python 3
7. Python libraries: Biopython https://biopython.org/

### Installation
Through this guide we assume the home directory is ??? and everything is installed there
Where to install,
How to install (external apps)


### Usage
This extension contains two main programs, optical map generator (???.py) and read colorer (???.py).

Optical map generator is used with 
```
python ???.py a b c d
```
a explanation
b explanation
c explanation
d explanation

Read colorer is used with 
```
python ???.py a b c d
```
a explanation
b explanation
c explanation
d explanation

### Example run



