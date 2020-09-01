#!/bin/bash

source activate ../envs/bowtie2

python mnutils.py -i ../mapping/WTV_default.sort.bam --scaffold scaffold51_166 --dump
# python mnutils.py -i ../mapping/WTV_default.sort.bam -o all --dump
