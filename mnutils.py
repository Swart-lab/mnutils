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
                    help="Calculate global phaseogram? (slow)")
parser.add_argument("--gff", type=str,
                    help="GFF3 file with features of interest")
parser.add_argument("--feature", type=str, default="five_prime_UTR",
                    help="Type of feature to produce phaseogram for")
args = parser.parse_args()

# Logging config
logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO)



def get_midpoint(segment):
    """
    Get presumed nucleosome dyad midpoint from digest fragment. 

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
    # midpoint = round(midpoint)
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
    filename : str
        Path to file to write PNG file
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
    """
    # Sort input 
    indict_sort = sorted(indict.items(), key=lambda item: int(item[0]))
    # xx = list(indict.keys())
    # yy = [indict[k] for k in xx]
    # xx = [int(i) for i in xx]
    # df['phase'] = [int(i) for i in df.index]
    xx = [int(i) for i,j in indict_sort]
    yy = [int(j) for i,j in indict_sort]

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
    plt.plot(xx,yy)
    plt.savefig(filename)


def feature_starts_from_gff3(filename, target_feature="five_prime_UTR"):
    """Parse start positions and orientations of feature type of interest from
    GFF3 file

    Parameters
    ----------
    filename : str
        Path to GFF3 file
    target_feature : str
        Type of feature to extract

    Returns
    -------
    defaultdict
        Dict of start positions keyed by scaffold -> orientation. Assumes that 
        features have orientation!
    """
    out = defaultdict(lambda: defaultdict(list))
    with open(filname, "r") as fh:
        for line in fh:
            line.rstrip()
            lsplit = line.split("\t")
            if len(lsplit) == 9: # Check for correctly formatted GFF3
                scaffold = lsplit[0]
                feature = lsplit[2]
                start = lsplit[3]
                stop = lsplit[4]
                orientation = lsplit[6]
                if feature == target_feature:
                    if orientation == "+":
                        out[scaffold]['+'].append(int(start))
                    elif orientation == "-":
                        out[scaffold]['-'].append(int(stop))
                        # We can do this because GFF coordinates are both 
                        # inclusive
                    else:
                        logging.warning(f"Feature {target_feature} at position {scaffold} {start} {stop} has invalid orientation {orientation}")
    return(out)


# -----------------------------------------------------------------------------

class Posmap(object):
    def __init__(self):
        """Construct new Posmap object
        
        _midpoints : defaultdict
            dict of nucleosome midpoints parsed from MNase-Seq mapping
            keyed by scaffold -> pos -> list of fragment lengths
        _tlen_histogram : defaultdict
            Histogram of fragment lengths
        """
        self._midpoints = defaultdict(lambda: defaultdict(list))
        self._tlen_histogram = defaultdict(int)


    def read_bam(self, bamfile, 
            *, scaffold=None, min_tlen=126, max_tlen=166):
        """Read BAM file and process it into nucleosome midpoint position map
        and fragment length histogram

        Parameters
        ----------
        bamfile : str
            Path to indexed BAM file to open with pysam
        scaffold : str
            Optional - parse this scaffold only
        min_tlen : int
            Lower length limit to consider a nucleosome fragment
        max_tlen : int
            Upper length limit to consider a nucleosome fragment
        """
        # Open BAM files
        logging.info(f"Opening BAM file {bamfile}")
        bamfh = pysam.AlignmentFile(bamfile, "rb")
        if scaffold:
            logging.info(f"Fetching only scaffold {scaffold}")
            bamiter = bamfh.fetch(scaffold)
        else:
            logging.info(f"Fetching all scaffolds")
            bamiter = bamfh.fetch()
        logging.info(f"Nucleosome length limits {min_tlen} and {max_tlen}")
        for read in bamiter:
            if (not read.is_secondary and
                    not read.is_duplicate and
                    read.get_tag("NM") <= 1):
                # TODO: Mismatches should actually be summed for each read PAIR,
                # not for each read segment as is currently done here
                if (read.template_length > 0): 
                    # To avoid double counting of read pairs, only consider 
                    # template lengths > 0, because template_length is negative
                    # for read2 in segment, and is 0 for unpaired
                    self._tlen_histogram[read.template_length] += 1
                if (read.template_length >= args.min_tlen and
                        read.template_length <= args.max_tlen and 
                        read.is_read1):
                    # Record midpoints only for reads within fragment length 
                    # criteria, and only for read1 in segment to avoid double
                    # counting
                    ref, pos, tlen = get_midpoint(read)
                    self._midpoints[ref][pos].append(int(tlen))
        bamfh.close

    def high_acc_frag_pc(self, nucl=147, window=2):
        """Find percentage of fragments with lengths close to true 
        nucleosome fragment length

        Parameters
        ----------
        nucl : int
            True length of nucleosome, in bp
        window : int
            Window (one-sided) of lengths to consider as "high accuracy"

        Returns
        -------
        float
            Percentage of fragments within high-accuracy length range
        """
        highacc = [int(self._tlen_histogram[k]) 
                    for k in range(nucl-window, nucl+1+window) 
                    if self._tlen_histogram[k]]
        tlen_total = sum([int(v) for v in self._tlen_histogram.values()])
        highacc_pc = 100 * sum(highacc) / tlen_total
        return(highacc_pc)


# MAIN # -----------------------------------------------------------------------

posmap = Posmap()

posmap.read_bam(args.input, 
        scaffold=args.scaffold, 
        min_tlen=args.min_tlen,
        max_tlen=args.max_tlen)

highacc_pc = posmap.high_acc_frag_pc()
logging.info(f"High accuracy fragments constitute {round(highacc_pc,2)}% of fragments")
if highacc_pc < 20:
    logging.info("Percentage is lower than 20%, high-accuracy method may not work well")

logging.info("Writing output files")
with open(f"{args.output}.tlen_hist.json","w") as fh:
    json.dump(posmap._tlen_histogram, fh, indent=4)
dict2plot_x_keys(posmap._tlen_histogram, title="Template length histogram", 
        xlabel="length", ylabel="counts", 
        xlim=(0,500), # hard-windowed to these limits for Illumina
        filename=f"{args.output}.tlen_hist.png")
with open(f"{args.output}.midpoints.json","w") as fh:
    json.dump(posmap._midpoints, fh, indent=4)

if args.phaseogram:
    logging.info("Calculating phaseogram")
    # TODO: Implement phaseogram calculation with midpoints

logging.info("Done")
