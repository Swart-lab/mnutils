#!/usr/bin/env python

# import re
import pysam
from collections import defaultdict
import argparse
import json
import logging
# import pandas as pd
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input BAM file")
parser.add_argument("-o", "--output", type=str,
                    help="Output prefix", default="test")
parser.add_argument("--scaffold", type=str, help="Scaffold to filter")
parser.add_argument("--min_tlen", type=int, default=126,
                    help="Minimum template insert length for read pair")
parser.add_argument("--max_tlen", type=int, default=166,
                    help="Maximum template insert length for read pair")
parser.add_argument("--max_mismatch", type=int, default=1,
                    help="Maximum mismatch allowed per alignment")
parser.add_argument("--phaseogram", 
                    help="Calculate phaseogram? (slow)")
args = parser.parse_args()

# Logging config
logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO)


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
        return(segment.reference_name, segment.reference_end, "rev")
    else:
        return(segment.reference_name, segment.reference_start, "fwd")


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
    with open(filename, "w") as fh:
        for tup in inlist:
            outstr = "\t".join([str(i) for i in tup])
            fh.write(outstr)
            fh.write("\n")


def fragstarts_dict_to_phaseogram(indict):
    """
    Convert dictionary of fragment start positions to phaseogram

    Parameters
    ----------
    indict : defaultdict
        Dict keyed by scaffolds -> pos -> frag start counts

    Returns
    -------
    list
        List of tuples (displacement, count)
    """

    # TODO: Normalize per scaffold
    # TODO: Average counts over total read start sites
    # TODO: Count frag starts upstream of position instead of downstream

    phaseogram = defaultdict(int)

    for scaffold in indict:
        for pos in indict[scaffold]:
            for i in range(1000):
                jumppos = i + 1 + int(pos)
                if jumppos in indict[scaffold].keys():
                    phaseogram[i+1] += indict[scaffold][jumppos]

    return(phaseogram)


def dict2plot_x_keys(indict, filename, 
        *, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None):
    """
    Plot phaseogram to PNG file

    Parameters
    ----------
    indict : dict
        dict containing parameters to plot. keys should be convertible to int.
        keys will be used as x-axis, values will be used as y-axis.
    title : str
        Title for the plot
    xlabel : str
        Label for x-axis
    ylabel : str
        Label for y-axis
    xlim : tuple
        tuple of limits (start, stop) for x-axis
    ylim : tuple
        tuple of limits (start, stop) for y-axis
    filename : str
        Path to file to write PNG file
    """
    # df = pd.DataFrame.from_dict(phaseogram, orient="index", columns=['count'])
    xx = list(indict.keys())
    yy = [indict[k] for k in xx]
    # df['phase'] = [int(i) for i in df.index]
    xx = [int(i) for i in xx]

    plt.figure()
    if title:
        plt.title(title)
    if xlabel:
        plt.xlabel(xlabel)
    if ylabel:
        plt.ylabel(ylabel)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    # plt.scatter(df['phase'], df['count'])
    plt.scatter(xx,yy)
    plt.savefig(filename)


# -----------------------------------------------------------------------------

# Open BAM files
logging.info(f"Opening SAM file {args.input}")
samfile = pysam.AlignmentFile(args.input, "rb")
# samout = pysam.AlignmentFile(args.output,"wb",template=samfile)

if args.scaffold:
    logging.info(f"Fetching only scaffold {args.scaffold}")
    samiter = samfile.fetch(args.scaffold)
else:
    logging.info(f"Fetching all scaffolds")
    samiter = samfile.fetch()

# Dict of scaffold -> pos -> frag start counts
fragstarts_fwd = defaultdict(lambda: defaultdict(int))
fragstarts_rev = defaultdict(lambda: defaultdict(int))

# List of tuples of mipdoints
midpoints = []

# Length histogram
tlen_histogram = defaultdict(int)

# for each aligned segment, report only those segments with template lengths
# within min/max limits
logging.info("Iterating through all reads in SAM file")
logging.info(f"Nucleosome length limits {args.min_tlen} and {args.max_tlen}")
for read in samiter:
    if (not read.is_secondary and
            not read.is_duplicate and
            read.get_tag("NM") <= 1):
        # TODO: Filter out reads with more than one mismatch, with NM tag
        # Actually mismatches should be summed for each read PAIR, not just
        # each read

        # Count number of fragment starts per base position per scaffold
        ref, pos, ori = get_frag_start(read)
        if ori == "fwd":
            fragstarts_fwd[ref][pos] += 1
        elif ori == "rev":
            fragstarts_rev[ref][pos] += 1

        if (read.template_length > 0): 
            # Avoid double counting, template_length is negative for rev reads
            # and is zero for unpaired
            tlen_histogram[read.template_length] += 1

        if (read.template_length >= args.min_tlen and
                read.template_length <= args.max_tlen):
            # Record midpoints only for reads within fragment length criteria
            # samout.write(read)
            midpoints.append(get_midpoint(read))

samfile.close

# listtuples_to_tsv(fragstarts,"test_fragstarts.tsv")
logging.info("Writing output files")
with open(f"{args.output}.fragstarts_fwd.json", "w") as fh:
    json.dump(fragstarts_fwd, fh, indent=4)
with open(f"{args.output}.fragstarts_rev.json", "w") as fh:
    json.dump(fragstarts_rev, fh, indent=4)
with open(f"{args.output}.tlen_hist.json","w") as fh:
    json.dump(tlen_histogram, fh, indent=4)
dict2plot_x_keys(tlen_histogram, title="Template length histogram", 
        xlabel="length", ylabel="counts", 
        xlim=(0,500),
        filename=f"{args.output}.tlen_hist.png")
listtuples_to_tsv(midpoints, f"{args.output}.midpoints.tsv")

if args.phaseogram:
    logging.info("Calculating phaseogram")
    phaseogram_fwd = fragstarts_dict_to_phaseogram(fragstarts_fwd)
    # phaseogram_rev = fragstarts_dict_to_phaseogram(fragstarts_rev)
    with open(f"{args.output}.phaseogram_fwd.json", "w") as fh:
        json.dump(phaseogram_fwd, fh, indent=4)
    # with open("test_phaseogram_rev.json", "w") as fh:
    #     json.dump(phaseogram_rev, fh, indent=4)

logging.info("Done")
