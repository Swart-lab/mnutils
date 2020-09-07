#!/usr/bin/env python3

import pysam
import logging
from mnutils import Shared
from collections import defaultdict
from random import sample


class Posmap(object):
    def __init__(self):
        """Construct new Posmap object
        
        _posdict : defaultdict
            dict of nucleosome positions parsed from MNase-Seq mapping
            keyed by scaffold -> pos -> list of fragment lengths
        _positionmap : dict
            dict of nucleosome raw position coverage
            keyed by scaffold -> pos -> no. nucleosome fragments
        _positionmap_smooth : dict
            dict of smoothed position coverage
            keyed by scaffold -> list of position coverage values
        _tlen_histogram : defaultdict
            Histogram of fragment lengths
        """
        self._posdict = defaultdict(lambda: defaultdict(list))
        self._positionmap = dict()
        self._positionmap_list = dict()
        self._positionmap_smooth = dict()
        self._tlen_histogram = defaultdict(int)


    def read_bam(self, bamfile,
            *, scaffold=None, min_tlen=126, max_tlen=166):
        """Read BAM file and process it into nucleosome position map and 
        fragment length histogram

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
                if (read.template_length >= min_tlen and
                        read.template_length <= max_tlen and 
                        read.is_read1):
                    # Record midpoints only for reads within fragment length 
                    # criteria, and only for read1 in segment to avoid double
                    # counting
                    ref, pos, tlen = Shared.get_midpoint(read)
                    self._posdict[ref][pos].append(int(tlen))
        bamfh.close
        logging.info("Calculating nucleosome position map")
        for ref in self._posdict:
            self._positionmap[ref] = {int(i): len(j) for i, j in self._posdict[ref].items()}


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


    def smooth_gaussian_positionmap(self, windowsize=50, bandwidth=10):
        """Gaussian smooth position map

        Parameters
        ----------
        windowsize : int
            Size of sliding window for smoother
        bandwidth : int
            Radius of Gaussian smoother

        Output to self._positionmap_smooth is a dict of lists, keyed by 
        scaffold. Each list contains array of smoothed values, with interval 1
        """
        self._positionmap_smooth = {}
        for scaffold in self._positionmap:
            yy_raw = [] # Raw unsmoothed values
            max_x = max([int(i) for i in self._positionmap[scaffold]])
            # Convert position map dict to array, filling in zeroes
            for i in range(1,max_x+1): # Start from 1 because 1 -based numbering
                if i in self._positionmap[scaffold]:
                    yy_raw.append(self._positionmap[scaffold][i])
                else:
                    yy_raw.append(0)
            yy_smooth, x_smooth_start = Shared.smooth_gaussian_integer_intervals(yy_raw, 1,1, 50, 10)
            # Prepend zeroes to smoothed curve
            yy_smooth = [0.0] * int(x_smooth_start) + yy_smooth
            self._positionmap_smooth[scaffold] = yy_smooth


    def smooth_gaussian_positionmap_to_wig(self, filename):
        """Write Wiggle Track Format (WIG) file of smoothed position map

        Parameters
        ----------
        filename : str
            Path to write WIG file
        """
        if len(self._positionmap_smooth) == 0:
            warnings.warn("Smoothed position map was not calculated")
        with open(filename, "w") as fh:
            for scaffold in self._positionmap_smooth:
                # WIG files use 1-based numbering
                # TODO: Check for off-by-one errors...
                fh.write(f"fixedStep chrom={scaffold} start=1 step=1\n")
                for i in self._positionmap_smooth[scaffold]:
                    fh.write(f"{i}\n")


    def global_phaseogram(self, window=1000, subsample=1000):
        """Calculate phaseogram of nucleosome positions up and downstream
        from other nucleosome.

        Parameters
        ----------
        window : int
            Window (symmetrical, up and downstream) in bp
        subsample : int
            Randomly sample this number of positions per scaffold.
            If the number of positions is less, then use all positions.
        """
        phaseogram = defaultdict(int)
        phaseogram[0] = 0

        for scaffold in self._positionmap:
            if len(self._positionmap[scaffold]) > subsample:
                poss = sample(list(self._positionmap[scaffold]), subsample)
            else:
                poss = self._positionmap[scaffold]
            for pos in poss:
                left = int(pos - window)
                right = int(pos + window)
                width = right - left
                for j in range(2 * window):
                    jumppos = left + j
                    if jumppos in self._positionmap[scaffold]:
                        phaseogram[j - window] += self._positionmap[scaffold][jumppos]
        return(phaseogram)


    def feature_phaseogram(self, gff3_file, 
            target_feature="five_prime_UTR", window=1000):
        """Calculate phaseogram of nucleosome positions around start positions
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
        feature_starts = Shared.feature_starts_from_gff3(gff3_file, target_feature)
        phaseogram = defaultdict(int)
        for scaffold in feature_starts:
            if scaffold in self._positionmap:
                for orientation in feature_starts[scaffold]:
                    for fwd_start in feature_starts[scaffold][orientation]:
                        left = fwd_start - window
                        right = fwd_start + window
                        width = right - left
                        for j in range(width):
                            if orientation == '+':
                                jumppos = left + j
                            elif orientation == '-':
                                jumppos = right - j
                            if jumppos in self._positionmap[scaffold]:
                                phaseogram[j - window] += self._positionmap[scaffold][jumppos]
            else:
                logging.warning(f"Scaffold {scaffold} in GFF file {gff3_file} but not in mapping")
        return(phaseogram)


