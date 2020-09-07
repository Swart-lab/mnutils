#!/usr/bin/env python

import pysam
from collections import defaultdict
import argparse
import json
import logging
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
parser.add_argument("--locmap", action="store_true",
                    help="Calculate localization measure sensu Zhang et al.")
parser.add_argument("--phaseogram", action="store_true",
                    help="Calculate global phaseogram? (slow)")
parser.add_argument("--gff", type=str,
                    help="GFF3 file with features of interest")
parser.add_argument("--feature", type=str, default="five_prime_UTR",
                    help="Type of feature to produce phaseogram for")
args = parser.parse_args()

# Logging config
logging.basicConfig(format='[%(asctime)s] %(message)s', level=logging.INFO)

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
    Shared.dict2plot_x_keys(phaseogram_feature, title=f"Phaseogram around {args.feature} starts",
            xlabel="Position vs. feature start (bp)",
            ylabel="Nucleosome position counts",
            xlim=(-1000,1000),
            filename=f"{args.output}.phaseogram.{args.feature}.png")
    if args.dump:
        with open(f"{args.output}.phaseogram.{args.feature}.json","w") as fh:
            json.dump(phaseogram_feature, fh, indent=4)

logging.info("Smoothing position map")
posmap.smooth_gaussian_positionmap(windowsize=50, bandwidth=10)

if args.locmap:
    logging.info("Calculating localization score sensu Zhang et al.")
    posmap.calculate_locmap()
    posmap.smooth_gaussian_locmap()

logging.info("Writing output files")
Shared.dict2plot_x_keys(posmap._tlen_histogram, title="Template length histogram", 
        xlabel="Length (bp)", ylabel="Counts", 
        xlim=(0,500), # hard-windowed to these limits for Illumina
        filename=f"{args.output}.tlen_hist.png")
posmap.write_wig(f"{args.output}.posmap_raw.wig","raw")
posmap.write_wig(f"{args.output}.posmap_smooth.wig","smooth")
if args.locmap:
    posmap.write_wig(f"{args.output}.locmap_raw.wig","raw_locmap")
    posmap.write_wig(f"{args.output}.locmap_smooth.wig","smooth_locmap")

if args.dump:
    logging.info("Dumping internal data to JSON files")
    with open(f"{args.output}.tlen_hist.json","w") as fh:
        json.dump(posmap._tlen_histogram, fh, indent=4)
    with open(f"{args.output}.posdict.json","w") as fh:
        json.dump(posmap._posdict, fh, indent=4)
    with open(f"{args.output}.posmap.json","w") as fh:
        json.dump(posmap._positionmap, fh, indent=4)
    with open(f"{args.output}.posmap_raw.json", "w") as fh:
        json.dump(posmap._positionmap_list, fh, indent=4)
    with open(f"{args.output}.posmap_smooth.json", "w") as fh:
        json.dump(posmap._positionmap_smooth, fh, indent=4)
    if args.locmap:
        with open(f"{args.output}.locmap_raw.json", "w") as fh:
            json.dump(posmap._locmap_list, fh, indent=4)
        with open(f"{args.output}.locmap_smooth.json", "w") as fh:
            json.dump(posmap._locmap_smooth, fh, indent=4)

if args.phaseogram:
    logging.info("Calculating global phaseogram with subsample of 5000 per scaffold")
    phaseogram_global = posmap.global_phaseogram(subsample=5000)
    Shared.dict2plot_x_keys(phaseogram_global, title=f"Global phaseogram for nucleosome positions",
            xlabel="Position vs. nucleosome (bp)",
            ylabel="Nucleosome position counts",
            xlim=(-1000,1000),
            filename=f"{args.output}.phaseogram_global.png")
    if args.dump:
        with open(f"{args.output}.phaseogram_global.json","w") as fh:
            json.dump(phaseogram_global, fh, indent=4)

logging.info("Done")
