#!/usr/bin/env python

from math import e
import warnings
import logging
import matplotlib.pyplot as plt
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
    midpoint = int(midpoint + 0.5) 
    # builtin round(midpoint) doesn't work properly for some reason, will 
    # always yield an even number, e.g. round(1.5) -> 2 and round(4.5) -> 4
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


def dict_to_wig(indict, filename):
    """Write Wiggle Track Format (WIG) file from dict of lists

    Parameters
    ----------
    indict : dict
        dict of lists, where keys are scaffold names
    filename : str
        Path to write WIG file
    """
    with open(filename, "w") as fh:
        for scaffold in indict:
            # WIG files use 1-based numbering
            # TODO: Check for off-by-one errors...
            fh.write(f"fixedStep chrom={scaffold} start=1 step=1\n")
            for i in indict[scaffold]:
                fh.write(f"{i}\n")


def localization_measure(yy_raw, int_rad=20, ext_rad=80):
    """Calculate the localization measure of Zhang et al. from raw nucleosome
    position map.

    Parameters
    ----------
    yy_raw : list
        Raw nucleosome position map, i.e. counts of nucleosomes per position for
        each position, at 1 bp intervals.
    int_rad : int
        Internal window radius. Note that this is a radius so the internal 
        window width is 2 * int_rad
    ext_rad : int
        External window radius. Note that this is a radius, like int_rad

    Returns
    -------
    list
        Localization measure values per position, at 1 bp intervals
    """
    out = [0.0] * ext_rad # Initialize output with zeroes for the left margin
    for center in range(ext_rad, len(yy_raw) - ext_rad):
        int_window = yy_raw[center - int_rad : center + int_rad]
        ext_window = yy_raw[center - ext_rad : center + ext_rad]
        if sum(ext_window) == 0:
            loc = 0
        else:
            loc = float(sum(int_window) / sum(ext_window))
        out.append(loc)
    out.extend([0.0] * ext_rad) # Pad right margin with zeroes
    return(out)


def find_peaks_iter(yy, mask_rad=70, min_peak=0.1):
    """Iteratively find peaks and mask surrounding radius until no more peaks
    found

    Parameters
    ----------
    yy : list
        list of values to search for peaks
    mask_rad: int
        Radius around found peaks to be masked; peaks must be at least this
        distance apart. NB this is a radius so window is 2*mask_rad
    min_peak : float
        Minimum (exclusive) value for peak to be assigned. Algorithm stops when 
        values in list are not higher than this.
    """
    yy_copy = yy.copy() # Make copy of the original list that we can copy
    out=[] # record peak positions
    out2=[] # for peaks whose radii overlap with existing peaks
    while(max(yy_copy) > min_peak):
        # Get highest point which we call a peak
        peak_idx = yy.index(max(yy))
        if min(yy_copy[peak_idx - mask_rad : peak_idx + mask_rad]) > 0:
            # If peak is not within radius of another peak, append to list
            out.append(peak_idx)
            # Mark radius around peak as -1
            yy_copy[peak_idx - mask_rad : peak_idx + mask_rad] = [-1] * (2*mask_rad)
        else:
            # Look for peaks that are within radius of another peak
            # Append to a separate list
            out2.append(peak_idx)
            yy_copy[peak_idx - mask_rad : peak_idx + mask_rad] = [-2] * (2*mask_rad)
    return(out, out2)


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



