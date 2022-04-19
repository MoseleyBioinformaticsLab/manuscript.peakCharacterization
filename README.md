This is a code repository to reproduce the following manuscript:

Scan-Centric, Frequency-Based Method for Characterizing Peaks from Direct Injection Fourier transform Mass Spectrometry Experiments

Robert M Flight, Joshua M Mitchell, Hunter N.B. Moseley

doi: https://doi.org/10.1101/2022.04.14.488423

```
git clone https://github.com/MoseleyBioinformaticsLab/manuscript.peakCharacterization.git

cd manuscript.peakCharacterization
```

We also need the data files and the _targets directory.

```
wget https://zenodo.org/record/6468443/files/scancentric_manuscript_targets.zip
unzip scancentric_manuscript_targets.zip

wget https://zenodo.org/record/6462315/files/scancentric_manuscript_data.zip
unzip scancentric_manuscript_data.zip
```

Start R within the manuscript directory:

```r
# depending on how it starts up, you may
# need to install renv first
install.packages("renv")
renv::install("BiocManager@1.30.16")
renv::install("broom@0.7.11")
renv::install("tictoc@1.0.1")
renv::install("bioc::xcms@3.16.1")
renv::restore()
```

After restoring the data, _targets, and R packages, we are ready to see if everything has been maintained.

If you look closely, you will notice that a lot of the data in the _targets.R file and other files is not listed here.
That is because it is 5.1 GB after zipping, so there was no way we were putting that in the git repo.

There will be instructions on how to clone, setup and get this running, just not today.
I need a break after getting the manuscript done.
