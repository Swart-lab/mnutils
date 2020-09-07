Nucleosome positioning utilities
--------------------------------

Software utility for QC and evaluation of nucleosome positioning data from
MNase-Seq or DNase-Seq experiments.

# Input data

 * BAM/SAM mapping of MNase-Seq reads to reference genome. Original sequencing
   data must have been paired-end sequencing.
 * Fasta file of reference genome sequence
 * GFF file of genomic features corresponding to sequence file

# Output data

 - [x] DNA fragment length histogram, evaluates effectiveness of digestion
 - [ ] Nucleosome occupancy map (simply coverage pileup for specific fragment range)
 - [x] Nucleosome absolute positioning map (using nucleosome _centers_ rather than 
   read starts)
   - [ ] With the high-accuracy procedure of Cole (2012)
   - [x] With broader fragment range
   - [x] Smoothed positioning map
   - [ ] Occupancy map (smoothing with uniform kernel and width 73)
 - [ ] Nucleosome conditional positioning map (like Kaplan 2010)
 - [x] Nucleosome global phaseogram (using nucleosome _centers_ vs read starts)
 - [x] Nucleosome localization measure map (using procedure of Zhang et al.)
 - [ ] Nucleosome calls (using iterative procedure of Kaplan et al 2010)
 - [x] Nucleosome feature-specific phaseogram (e.g. relative to TSS)
 - [ ] Well positioned nucleosome array calls (using FFT procedure)
