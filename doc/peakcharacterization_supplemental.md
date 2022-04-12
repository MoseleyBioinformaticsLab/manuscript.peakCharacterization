---
title: "Supplemental Materials for: Scan-Centric, Frequency-Based Method for Characterizing Peaks from Direct Injection Fourier transform Mass Spectrometry Experiments"
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
date: "2022-04-12 11:38:43"
output: 
  word_document:
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

## Bruker SolariX ICR Frequency Conversion

In contrast to the m/z to frequency conversion equations used for Thermo-Fisher Fusion and other Orbitrap instruments, the one we have used for Bruker SolariX ICR instruments is simpler:

$$frequency = a + x * \frac{1}{mz} + y * \frac{1}{\sqrt{mz}}$$
