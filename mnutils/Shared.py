#!/usr/bin/env python

from math import e
import warnings
import logging
from collections import defaultdict


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
    # Take absolute because template_length can be negative for reverse reads
    tlen = abs(segment.template_length)
    if segment.is_reverse:
        midpoint = segment.reference_end - (tlen/2)
    else:
        midpoint = segment.reference_start + (tlen/2)
    # midpoint = round(midpoint)
    return(ref, midpoint, tlen)


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


def gaussian_kernel_coeffs(bandwidth, xx, x_null):
    """Precompute coefficients for Gaussian smoothing kernel
    
    Parameters
    ----------
    bandwidth : float
        Smoothing bandwidth
    xx : list
        List of floats of X-values for which to compute the coefficients.
    x_null : float
        Center X-value

    Return
    ------
    list
        List of coefficients K corresponding to the input X-values
    """

    out = [pow(e, -0.5 * (x - x_null)**2 / bandwidth**2) for x in xx]

    return(out)


# TODO: Implement uniform kernel and generalize smoothing function


def smooth_gaussian_integer_intervals(yy_raw, x_start, x_interval, x_window, bandwidth):
    """Perform Gaussian smoothing on a list of raw Y-values where the 
    X-values are of constant integer intervals

    Parameters
    ----------
    yy_raw : list
        List of floats, raw Y-values
    x_start : int
        X-value corresponding to the initial Y-value
    x_interval : int
        Interval between Y-values 
    x_window : int
        Window (kernel) width, should be an even number and about >3-fold
        the bandwidth
    bandwidth : float
        Bandwidth for the Gaussian kernel function

    Returns
    -------
    list, int
        List of smoothed Y-values, and starting X-value
    """
    if x_interval > x_window:
        warnings.warn("Step interval cannot be larger than window for smoothing")
    if x_window < bandwidth:
        warnings.warn("Bandwidth should be smaller than kernel window")
    if x_window % 2 != 0:
        warnings.warn("Kernel window should be an even integer")
    # Number of values per window
    num_intervals_half_window = int(x_window / (2 * x_interval))
    N = int((x_window / x_interval ) + 1)
    # Get window of X-values centered on 0
    xx = [0 - x_window/2 + k * x_interval for k in range(N)]
    # Caclulate kernel coefficients for those X-values
    coeff_K = gaussian_kernel_coeffs(bandwidth, xx, 0)
    norm_K = [k / sum(coeff_K) for k in coeff_K]
    # Initialize output list
    yy_hat = [] # Smoothed Y-values

    for i in range(len(yy_raw)):
        x_null = x_start + i*x_interval
        # Start after half-window from beginning
        if i >= num_intervals_half_window and i < (len(yy_raw) - num_intervals_half_window):
            # Get window of Y-values
            yy_window = yy_raw[int(i-num_intervals_half_window) : int(i+num_intervals_half_window) + 1]
            # lengths of yy_window, xx, and coeff_K should be the same
            # Get weighted Y values by multiplying against normalized coefficients
            yy_adj = [y * k for y, k in zip(yy_window, norm_K)]
            yy_hat.append(sum(yy_adj))

    return(yy_hat, x_start + x_window/2)


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
