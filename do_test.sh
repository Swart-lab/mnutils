#!/bin/bash

source activate ../envs/bowtie2

python mnutils.py -i ../mapping/WTV_default.sort.bam --scaffold scaffold51_166
