---
title: "Scan-Centric, Frequency-Based Method for Characterizing Peaks from Direct Injection Fourier transform Mass Spectrometry Experiments"
author:
  - Robert M Flight:
      email: robert.flight@uky.edu
      institute: [markey, biochem, rcsirm]
  - Joshua M Mitchell:
      email: jmmi243@uky.edu
      institute: [markey, biochem, rcsirm, ibi]
  - Hunter NB Moseley:
      email: hunter.moseley@uky.edu
      correspondence: true
      institute: [markey, biochem, rcsirm, ibi, tox]
institute:
  - markey: Markey Cancer Center, University of Kentucky, Lexington, KY 40536, USA
  - biochem: Department of Molecular & Cellular Biochemistry, University of Kentucky, Lexington, KY 40536, USA
  - rcsirm: Resource Center for Stable Isotope Resolved Metabolomics, University of Kentucky, Lexington, KY 40536, USA
  - ibi: Institute for Biomedical Informatics, University of Kentucky, Lexington, KY 40536, USA
  - tox: Department of Toxicology and Cancer Biology, University of Kentucky, Lexington, KY 40536, USA
date: "2022-04-05 09:28:45"
output: 
  word_document:
    reference_docx: metabolites-template-v3.docx
    keep_md: true
    pandoc_args:
      - '--lua-filter=scholarly-metadata.lua'
      - '--lua-filter=author-info-blocks.lua'
bibliography: '/home/rmflight/Documents/manuscripts/in_progress/rmflight_peakCharacterization_new/doc/peakcharacterization.json'
csl: plos-computational-biology.csl
editor_options: 
  chunk_output_type: console
  markdown: 
    wrap: sentence
---





::: {custom-style="MDPI_1.7_abstract"}
**Abstract:** We introduce a novel, scan-centric method for characterizing peaks from direct injection multi-scan Fourier transform mass spectra of complex samples that utilizes frequency values derived directly from the spacing of raw m/z points in spectral scans.
Our peak characterization method utilizes intensity independent noise removal and normalization of scan level data to provide a much better fit of relative intensity to natural abundance probabilities for low abundance isotopologues that are not present in all of the acquired scans.
Moreover, our method calculates both peak- and scan-specific statistics incorporated within a series of quality control steps that are designed to robustly derive peak centers, intensities and intensity ratios with their scan-level variances.
These cross-scan characterized peaks are suitable for use in our previously published peak assignment methodology, Small Molecule Isotope Resolved Formula Enumeration (SMIRFE).
:::

::: {custom-style="MDPI_1.8_keywords"}
**Keywords**: Fourier-transform mass-spectrometry
:::

### Deleteme



::: {custom-style="MDPI_3.9_equation"}
$$
\begin{align}
frequency = intercept + x* \frac{1}{\sqrt{mz}} + y * \frac{1}{\sqrt[3]{mz}}&& \text{(1)}
\end{align}$$
:::


::: {custom-style="MDPI_3.9_equation"}
$$
\begin{align}
frequency = intercept + x* \frac{1}{\sqrt{mz}} + y * \frac{1}{\sqrt[3]{mz}} && \text{(1)}
\end{align}
$$
:::


::: {custom-style="MDPI_3.9_equation"}
*frequency* = a + x * *mz*^-1/2^ + y * *mz*^-1/3^
:::
::: {custom-style="MDPI_3.a_equation_number"}
(1)
:::

mz = a + x * *frequency*^-1^ + y * *frequency*^-2^ + z * *frequency*^-3^

ln(*intensity*) = a + x * *position* + y * *position*^2^

ln(*NAP*~P1~ / *NAP*~P2~) - ln(*Int*~P1~ / *Int*~P2~) \U2248 0
