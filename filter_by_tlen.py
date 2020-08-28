#!/usr/bin/env python

import re
import pysam
from collections import defaultdict
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i","--input",type=str,help="Input BAM file")
parser.add_argument("-o","--output",type=str,help="Output file",default="filter_out.bam")
parser.add_argument("--scaffold",type=str,help="Scaffold to filter")
parser.add_argument("--min_tlen",type=int,default=126,help="Minimum template insert length for read pair")
parser.add_argument("--max_tlen",type=int,default=166,help="Maximum template insert length for read pair")
parser.add_argument("--max_mismatch",type=int,default=1,help="Maximum mismatch allowed per alignment")
args = parser.parse_args()

def get_midpoint(segment):
    """
    Get presumed nucleosome dyad midpoint from digest fragment

    Parameters
    ----------
    segment : pysam.AlignedSegment
        Aligned read segment to parse. Should pass only read1 to avoid double
        counting of read pairs.

    Returns
    -------
    str, float, int
        tuple of reference name, midpoint pos, fragment length
    """
    ref = segment.reference_name
    if segment.is_reverse:
        midpoint = segment.reference_end - (segment.template_length/2)
    else:
        midpoint = segment.reference_start + (segment.template_length/2)
    return(ref, midpoint, segment.template_length)

def get_frag_start(segment):
    """
    Report fragment start position and orientation, for phaseogram

    Parameters
    ----------
    segment : pysam.AlignedSegment
        Aligned read segment to parse.

    Returns
    -------
    str, int, str
        tuple of reference name, fragment start pos, orientation
    """
    if segment.is_reverse:
        return(segment.reference_name,segment.reference_end,"rev")
    else:
        return(segment.reference_name, segment.reference_start,"fwd")

def listtuples_to_tsv(inlist, filename):
    """
    Write list of tuples to TSV file

    Parameters
    ----------
    inlist : list
        List of tuples to convert to TSV
    filename : str
        Path to file to write output
    """
    with open (filename, "w") as fh:
        for tup in inlist:
            outstr = "\t".join([str(i) for i in tup])
            fh.write(outstr)
            fh.write("\n")

#-------------------------------------------------------------------------------

# Open BAM files
samfile = pysam.AlignmentFile(args.input,"rb")
# samout = pysam.AlignmentFile(args.output,"wb",template=samfile)

if args.scaffold:
    samiter = samfile.fetch(args.scaffold)
else:
    samiter = samfile.fetch()

fragstarts = []
midpoints = []

# for each aligned segment, report only those segments with template lengths
# within min/max limits
for read in samiter:
    if not read.is_secondary and not read.is_duplicate and read.get_tag("NM") <= 1:
        # TODO: Filter out reads with more than one mismatch, with NM tag
        # Actually mismatches should be summed for each read PAIR, not just each read
        fragstarts.append(get_frag_start(read))
        if read.template_length >= args.min_tlen and read.template_length <= args.max_tlen:
            # samout.write(read)
            midpoints.append(get_midpoint(read))

samfile.close

listtuples_to_tsv(fragstarts,"test_fragstarts.tsv")
listtuples_to_tsv(midpoints,"test_midpoints.tsv")
