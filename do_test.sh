#!/bin/bash

source activate ../envs/bowtie2

# python mnutils.py -i ../mapping/WTV_default.sort.bam \
#	-o test \
# 	--scaffold scaffold51_166 \
# 	--gff scaffold51_166.gff3 \
# 	--dump

python mnutils.py -i ../mapping/WTV_default.sort.bam \
	-o all --dump \
	--gff ../ref/ptetraurelia_mac_51_annotation_v2.0.gff3
