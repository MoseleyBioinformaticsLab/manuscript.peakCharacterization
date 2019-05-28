# Outline

**Peak Characterization and Multi-Scan Averaging for FT-MS Data**

Note that this outline is more about figuring out the plots that need to be generated for this
manuscript. If we have an outline, then we can figure out the plots, and craft the
narrative around it.

## Abstract

## Introduction

* Multi-scan data is useful, b/c of failures in individual scans

* Common to simply average the scans together

* But is this correct?

* We will show it is not, and that an alternative is required

## Methods

* Raw data conversion

* Converting to frequency space
  * Models built using points that show consistent 0.5 difference in frequency space

* Removing scans
  * Based on outliers in their individual frequency model
* Conversion to offset or frequency space
  * tempting to use LOESS, but can't get a zero intercept
  * instead use points that show consistent 0.5 difference in frequency space

* Peak Characterization
  * up-down
  * linear model on the log-intensities

* Non-zero point densities
  * Window based on frequency
  * Sliding windows
  * Tiled windows
  * Remove < 99% for sliding windows (removes most of noise)
  * Reduce remaining
  * Find peaks in reduced
  * Count overlap with tiles
  * Split at zeros

* Normalize
  * Keep peaks > 0.7 log-ratio of max in a scan
  * Find most equidistant scan
  * Take median log-ratio to other scans
  * Normalize based on this ratio
  * Find peaks whose intensity is heavily correlated with scan number, remove them, re-normalize

## Results

### Demonstration of the Problem ...

* Need a result that shows the NAP ratio for something easily identifiable are way off
if you do simple averaging in XCalibur or in `xcms`

### Conversion to Frequency Space

* Plot showing how to convert to frequency space, with these sub-plots
  * Plot of M/Z point differences vs M/Z
  * Plot of frequency point differences vs frequency
  * Plot of frequency vs M/Z, and the analytical model fit
  * Plot of frequency intercept, and the outliers

### Removing Odd Scans

* Have outlier plot in the previous section
* May need to quantify across a data set how many get removed

### Peak Detection and Characterization

* Plot of raw spectrum showing low base-line
* Plot of center / intensity after doing peak characterization

### Noise Point Removal

* Histogram of non-zero counts for sliding regions

### Splitting Sub-Regions

* Plot of a multiple peak region
* Plot of counts in tiles

### Need for Noise Removal

* # of SMIRFE assignments with / without noise

### Need for Normalization

* # of SMIRFE assignments with / without normalization
* Show that normalization appears constant in a scan if use high log-ratio points

### Removing HPD Sites Using Frequency SD

* Find HPD sites in Xcalibur and `xcms` data, compare locations with final data before frequency sd peak
removal.

### Final Result

* NAP ratios now match expected, within some extremely tight tolerances, and M/Z matches rather well as well
* Also note that this removes the HPD regions observed previously
* Plot comparing RSD distributions / median values
  * With `noise` points left in
  * Without `noise` points left in
  * Without normalization
  * With single pass normalization using all peaks
  * With single pass normalization using 0.7 peaks
  * With double pass (removing scan correlated peaks) normalization
  * With high frequency SD peaks left in
  * With high frequency SD peaks removed

## Data Needed

* Xcalibur peak lists for all files
  * Josh very kindly generated these
* `xcms` averaged and our peak picking peak lists for all files
* combinations of characterized results needed
  * **group1**: noise kept, no normalization, no frequency sd filtering
  * **group2**: noise removed, no normalization, no frequency sd filtering
  * **group3**: noise removed, single pass normalization all peaks, no frequency sd filtering
  * **group4**: noise removed, single pass normalization 0.7 peaks, no frequency sd filtering
  * **group5**: noise removed, double pass normalization 0.7 peaks, no frequency sd filtering
  * **group6**: noise removed, double pass normalization 0.7 peaks, frequency sd filtering
