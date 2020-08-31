Nucleosome positioning utilities
--------------------------------

Software utility for QC and evaluation of nucleosome positioning data from
MNase-Seq or DNase-Seq experiments.

# Input data

 * BAM/SAM mapping of MNase-Seq reads to reference genome
 * Fasta file of reference genome sequence
 * GFF file of genomic features corresponding to sequence file

# Output data

 * DNA fragment length histogram, evaluates effectiveness of digestion
 * Nucleosome occupancy map (simply coverage pileup for specific fragment range)
 * Nucleosome absolute positioning map (using nucleosome _centers_ rather than 
   read starts)
   * With the high-accuracy procedure of Cole (2012)
   * With broader fragment range
 * Nucleosome conditional positioning map (like Kaplan 2010)
 * Nucleosome global phaseogram (using nucleosome _centers_ rather than read starts)
 * Nucleosome feature-specific phaseogram (e.g. relative to TSS)
