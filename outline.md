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

* Peak detection and peak characterization
  * up-down
  * linear model on the log-intensities

* Matching across scans
  * error model based on digital resolution
  * error model based on matching peaks
  * offset model

## Results

### Demonstration of the Problem ...

* Need a result that shows the NAP ratio for something easily identifiable are way off
if you do simple averaging in XCalibur (see #2)

### Peak Detection and Characterization

* Plot of raw spectrum showing low base-line
* Plot of center / intensity after doing peak characterization

### Noise Peak Removal

* Plot of log10 intensity

### Peak Matching Across Scans

* Plot of digital resolution model
* Plot of changes in error model
* Plot of changes in offsets

### Need for Noise Removal

* RSD is way off if they are included

### Need for Normalization

* RSD is still off
* Median noise levels
* Show that normalization is constant / not constant across M/Z (see #1)
* What do the NAP ratios look like w/out normalization?

### Final Result

* NAP ratios now match expected, within some extremely tight tolerances, and M/Z matches rather well as well
* Also note that this removes the HPD regions observed previously


