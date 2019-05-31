#!/bin/bash

# Target genome
#TARGET="ecoli"
#TARGET="scerevisiae"
TARGET="celegans"
#TARGET="arabidopsis"

# Provided data
REFERENCE="data/"$TARGET"-reference.fasta"
ERR_READS="data/"$TARGET"-reads.fasta"

# Correcting
COR_READS="data/"$TARGET"-reads-corrected.fasta"

# Other stuff

READ_READ_ALIGN_ZIP="data/"$TARGET"-reads-reads-alignment.paf.gz"
MINIASM_GFA="data/"$TARGET"-miniasm-output.gfa"
PRECOLOR_CONTIGS="data/"$TARGET"-precolor-contigs.fasta"

PRECOLOR_CONTIGS_OPT_MAP="data/"$TARGET"-precolor-contigs.valouev"
REFERENCE_OPT_MAP="data/"$TARGET"-reference.valouev"

CONTIG_REFERENCE_MAPPING="data/"$TARGET"-contig-reference-mapping-output.txt"

CONTIGS_READS_ALIGN="data/"$TARGET"-contigs-reads-alignment.paf"

COLORED_CONTIGS="data/"$TARGET"-precolor-contig-colors.cval"
COLORED_READS="data/"$TARGET"-read-colors.cf"
ADJ_COLORED_READS="data/"$TARGET"-read-colors-adjusted.cf"

KERMIT_GFA="data/"$TARGET"-kermit-output.gfa"
KERMIT_GFA_NP="data/"$TARGET"-kermit-np-output.gfa"

POSTCOLOR_CONTIGS="data/"$TARGET"-postcolor-contigs.fasta"
POSTCOLOR_CONTIGS_NP="data/"$TARGET"-postcolor-contigs-np.fasta"



###
### 1. Correct reads with CONSENT
###

#CONSENT/CONSENT-correct --in $ERR_READS --out $COR_READS --type PB

###
### 2. Run minimap with corrected reads 
###

#minimap2/minimap2 -x ava-pb -t8 $COR_READS $COR_READS | gzip -1 > $READ_READ_ALIGN_ZIP

###
### 3. Run miniasm to get contigs 
###

#miniasm/miniasm -f $COR_READS $READ_READ_ALIGN_ZIP > $MINIASM_GFA
#awk '$1 ~/S/ {print ">"$2"\n"$3}' $MINIASM_GFA | fold > $PRECOLOR_CONTIGS


###
### 4. Create optical maps of contigs
###

#python3 optical_map_generator.py $PRECOLOR_CONTIGS -o $PRECOLOR_CONTIGS_OPT_MAP


### 
### 5. Create optical map of reference
###

#python3 optical_map_generator.py $REFERENCE -o $REFERENCE_OPT_MAP


###
### 6. Run valouev to map contigs to reference 
###

valouev_optmap_alignment/fit/fit $REFERENCE_OPT_MAP $PRECOLOR_CONTIGS_OPT_MAP | tee $CONTIG_REFERENCE_MAPPING


###
### 7. Run minimap with corrected reads and contigs 
###

minimap2/minimap2 -t8 -x ava-pb $PRECOLOR_CONTIGS $COR_READS > $CONTIGS_READS_ALIGN


###
### 8.1 Color contigs
###

python3 contig_colorer.py $CONTIG_REFERENCE_MAPPING $REFERENCE_OPT_MAP $PRECOLOR_CONTIGS_OPT_MAP -c $COLORED_CONTIGS -t 40

###
### 8.2 Color reads based on contig colors and alignments
###

python3 read_colorer.py $CONTIGS_READS_ALIGN $COLORED_CONTIGS -r $COLORED_READS -e -s -l 0.5

###
### 9. Adjust colorings
###

python3 read_recolorer.py $COLORED_READS -a $ADJ_COLORED_READS

###
### 10. Run kermit
###

kermit/kermit -C $ADJ_COLORED_READS -f $COR_READS $READ_READ_ALIGN_ZIP > $KERMIT_GFA
kermit/kermit -C $ADJ_COLORED_READS -f $COR_READS $READ_READ_ALIGN_ZIP -P > $KERMIT_GFA_NP

awk '$1 ~/S/ {print ">"$2"\n"$3}' $KERMIT_GFA > $POSTCOLOR_CONTIGS
awk '$1 ~/S/ {print ">"$2"\n"$3}' $KERMIT_GFA_NP > $POSTCOLOR_CONTIGS_NP

python3 quast/quast-5.0.2/quast.py -r $REFERENCE $PRECOLOR_CONTIGS
python3 quast/quast-5.0.2/quast.py -r $REFERENCE $POSTCOLOR_CONTIGS
python3 quast/quast-5.0.2/quast.py -r $REFERENCE $POSTCOLOR_CONTIGS_NP




