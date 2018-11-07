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

* Removing scans
  * Probably based on TIC
* Conversion to offset or frequency space
  * tempting to use LOESS, but can't get a zero intercept
  * instead use points that show consistent 0.5 difference in frequency space

* Peak Characterization
  * up-down
  * linear model on the log-intensities

* Non-zero point densities
  * Window based on M/Z difference model
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

## Results

### Demonstration of the Problem ...

* Need a result that shows the NAP ratio for something easily identifiable are way off
if you do simple averaging in XCalibur (see #2)

### Removing Odd Scans

* From lung cancer data, maybe go through and count the original number of
scans acquired, and then how many scans on average get removed.

### Peak Detection and Characterization

* Plot of raw spectrum showing low base-line
* Plot of center / intensity after doing peak characterization

### Noise Point Removal

* Histogram of non-zero counts for sliding regions

### Splitting Sub-Regions

* Plot of a multiple peak region
* Plot of counts in tiles

### Need for Noise Removal

* RSD is way off if they are not removed
* # of SMIRFE assignments with / without noise

### Need for Normalization

* RSD is still off
* # of SMIRFE assignments with / without normalization
* Show that normalization appears constant in a scan if use high log-ratio points

### Final Result

* NAP ratios now match expected, within some extremely tight tolerances, and M/Z matches rather well as well
* Also note that this removes the HPD regions observed previously


