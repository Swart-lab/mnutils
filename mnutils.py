#!/usr/bin/env python

import pysam
from collections import defaultdict
import argparse
import json
import logging
import matplotlib.pyplot as plt
from mnutils import Shared, Posmap

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


# MAIN # -----------------------------------------------------------------------

posmap = Posmap.Posmap()

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
            ylabel="Nucleosome position counts",
            xlim=(-1000,1000),
            filename=f"{args.output}.phaseogram.{args.feature}.png")
    if args.dump:
        with open(f"{args.output}.phaseogram.{args.feature}.json","w") as fh:
            json.dump(phaseogram_feature, fh, indent=4)

logging.info("Smoothing position map")
posmap.smooth_gaussian_positionmap(windowsize=50, bandwidth=10)
posmap.smooth_gaussian_positionmap_to_wig(f"{args.output}.posmap_smooth.wig")

logging.info("Writing output files")
dict2plot_x_keys(posmap._tlen_histogram, title="Template length histogram", 
        xlabel="Length (bp)", ylabel="Counts", 
        xlim=(0,500), # hard-windowed to these limits for Illumina
        filename=f"{args.output}.tlen_hist.png")
if args.dump:
    logging.info("Dumping internal data to JSON files")
    with open(f"{args.output}.tlen_hist.json","w") as fh:
        json.dump(posmap._tlen_histogram, fh, indent=4)
    with open(f"{args.output}.posdict.json","w") as fh:
        json.dump(posmap._posdict, fh, indent=4)
    with open(f"{args.output}.posmap.json","w") as fh:
        json.dump(posmap._positionmap, fh, indent=4)
    with open(f"{args.output}.posmap_smooth.json", "w") as fh:
        json.dump(posmap._positionmap_smooth, fh, indent=4)

if args.phaseogram:
    logging.info("Calculating global phaseogram with subsample of 5000 per scaffold")
    phaseogram_global = posmap.global_phaseogram(subsample=5000)
    dict2plot_x_keys(phaseogram_global, title=f"Global phaseogram for nucleosome positions",
            xlabel="Position vs. nucleosome (bp)",
            ylabel="Nucleosome position counts",
            xlim=(-1000,1000),
            filename=f"{args.output}.phaseogram_global.png")
    if args.dump:
        with open(f"{args.output}.phaseogram_global.json","w") as fh:
            json.dump(phaseogram_global, fh, indent=4)

logging.info("Done")
