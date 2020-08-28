#!/bin/bash

source activate ../envs/bowtie2

python filter_by_tlen.py -i ../mapping/WTV.sort.bam --scaffold scaffold51_166
