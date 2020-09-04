#!/bin/bash

source activate ../envs/bowtie2

python mnutils.py -i ../mapping/WTV_default.sort.bam \
	-o testing_output/test2 \
	--scaffold scaffold51_166 \
	--gff scaffold51_166.gff3 \
	--phaseogram \
	--dump

# python mnutils.py -i ../mapping/WTV_default.sort.bam \
# 	-o testing_output/all2 \
# 	--gff ../ref/ptetraurelia_mac_51_annotation_v2.0.gff3 \
# 	--phaseogram \
# 	--dump 
