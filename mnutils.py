#!/usr/bin/env python

import pysam
from collections import defaultdict
from random import sample
import argparse
import json
import logging
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", type=str, help="Input BAM file")
parser.add_argument("-o", "--output", type=str,
                    help="Output prefix", default="test")
parser.add_argument("--dump", action="store_true",
                    help="Dump JSON files of internal objects")
parser.add_argument("--scaffold", type=str, help="Scaffold to filter")
parser.add_argument("--min_tlen", type=int, default=126,
                    help="Minimum template insert length for read pair")
parser.add_argument("--max_tlen", type=int, default=166,
                    help="Maximum template insert length for read pair")
parser.add_argument("--max_mismatch", type=int, default=1,
                    help="Maximum mismatch allowed per alignment")
parser.add_argument("--phaseogram", action="store_true",
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
        *, title=None, xlabel=None, ylabel=None, xlim=None, ylim=None,
        width=10, height=5):
    """
    Plot data of key-value pairs in a dict to PNG file

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
    width : int
        Width of plot (inches)
    height : int
        Height of plot (inches)
    """
    # Sort input 
    indict_sort = sorted(indict.items(), key=lambda item: float(item[0]))
    xx = [float(i) for i,j in indict_sort]
    yy = [float(j) for i,j in indict_sort]

    plt.figure(figsize=(width,height))
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
    with open(filename, "r") as fh:
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
        _positionmap : defaultdict
            dict of nucleosome midpoint position coverage
            keyed by scaffold -> pos -> no. nucleosome midpoints
        _tlen_histogram : defaultdict
            Histogram of fragment lengths
        """
        self._midpoints = defaultdict(lambda: defaultdict(list))
        self._positionmap = dict()
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
        logging.info("Calculating position map of nucleosome midpoints")
        for ref in self._midpoints:
            self._positionmap[ref] = {float(i): len(j) for i, j in self._midpoints[ref].items()}

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


    def global_phaseogram(self, window=1000, subsample=1000):
        """Calculate phaseogram of nucleosome midpoints up and downstream
        from nucleosome midpoints.

        Parameters
        ----------
        window : int
            Window (symmetrical, up and downstream) in bp
        subsample : int
            Randomly sample this number of positions per scaffold.
            If the number of positions is less, then use all positions.
        """
        phaseogram = defaultdict(int)

        for scaffold in self._positionmap:
            if len(self._positionmap[scaffold]) > subsample:
                poss = sample(list(self._positionmap[scaffold]), subsample)
            else:
                poss = self._positionmap[scaffold]
            for pos in poss:
                left = pos - window
                right = pos + window
                width = right - left
                for j in range(2 * 2 * window):
                    jumppos = left + j/2
                    if jumppos in self._positionmap[scaffold].keys():
                        phaseogram[j/2 - window] += self._positionmap[scaffold][jumppos]

        return(phaseogram)


    def feature_phaseogram(self, gff3_file, 
            target_feature="five_prime_UTR", window=1000):
        """Calculate phaseogram of nucleosome midpoints around start positions
        of features of interest.

        Parameters
        ----------
        gff3_file : str
            Path to GFF3 file containing annotations
        target_feature : str
            Type of feature of interest; should have +/- orientations
        window : int
            Width up- and down-stream of feature starts to make phaseogram
        """
        feature_starts = feature_starts_from_gff3(gff3_file, target_feature)
        phaseogram = defaultdict(int)
        for scaffold in feature_starts:
            if scaffold in self._positionmap:
                for orientation in feature_starts[scaffold]:
                    for fwd_start in feature_starts[scaffold][orientation]:
                        left = fwd_start - window
                        right = fwd_start + window
                        width = right - left
                        for j in range(2 * width):
                            jumppos = float()
                            if orientation == '+':
                                jumppos = left + j/2 # Account for half-integer positions
                            elif orientation == '-':
                                jumppos = right - j/2
                            if jumppos in self._positionmap[scaffold]:
                                phaseogram[j/2 - window] += self._positionmap[scaffold][jumppos]
            else:
                logging.warning(f"Scaffold {scaffold} in GFF file {gff3_file} but not in mapping")
        return(phaseogram)


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

if args.gff:
    logging.info(f"Parsing feature starts from GFF3 file {args.gff}")
    logging.info(f"Calculating phaseogram for {args.feature} features")
    phaseogram_feature = posmap.feature_phaseogram(args.gff, 
            target_feature=args.feature, 
            window=1000)
    dict2plot_x_keys(phaseogram_feature, title=f"Phaseogram around {args.feature} starts",
            xlabel="Position vs. feature start (bp)",
            ylabel="Nucleosome midpoint counts",
            xlim=(-1000,1000),
            filename=f"{args.output}.phaseogram.{args.feature}.png")
    if args.dump:
        with open(f"{args.output}.phaseogram.{args.feature}.json","w") as fh:
            json.dump(phaseogram_feature, fh, indent=4)

logging.info("Writing output files")
dict2plot_x_keys(posmap._tlen_histogram, title="Template length histogram", 
        xlabel="Length (bp)", ylabel="Counts", 
        xlim=(0,500), # hard-windowed to these limits for Illumina
        filename=f"{args.output}.tlen_hist.png")
if args.dump:
    logging.info("Dumping internal data to JSON files")
    with open(f"{args.output}.tlen_hist.json","w") as fh:
        json.dump(posmap._tlen_histogram, fh, indent=4)
    with open(f"{args.output}.midpoints.json","w") as fh:
        json.dump(posmap._midpoints, fh, indent=4)

if args.phaseogram:
    logging.info("Calculating global phaseogram with subsample of 5000 per scaffold")
    # TODO: Implement phaseogram calculation with midpoints
    phaseogram_global = posmap.global_phaseogram(subsample=5000)
    dict2plot_x_keys(phaseogram_global, title=f"Global phaseogram for nucleosome midpoints",
            xlabel="Position vs. nucleosome midpoint (bp)",
            ylabel="Nucleosome midpoint counts",
            xlim=(-1000,1000),
            filename=f"{args.output}.phaseogram_global.png")
    if args.dump:
        with open(f"{args.output}.phaseogram_global.json","w") as fh:
            json.dump(phaseogram_global, fh, indent=4)

logging.info("Done")
