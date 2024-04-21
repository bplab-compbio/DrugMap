# Cysteine Set Enrichment Analysis

This directory contains the code required to run CSEA (Cysteine Set Enrichment Analysis), as well as the initial repository of cysteine sets which were used in DrugMap.

<p align="center">
  <img src="https://github.com/bplab-compbio/DrugMap/blob/main/src/images/csea.png" width="1200" height="300">
</p>

At a glance, what does CSEA do?

1. You, the cysteine hunter, have a rarefied list of interesting cysteines which you want to situate in the broader context of the known biological universe.
2. Then, you ask whether these cysteines are common with cysteines with known biological roles (i.e. EGF-EGFR-PI3K signaling axis) **OR** cysteines which localize to a particular structural element (i.e. a pocket) **OR** with cysteines whose reactivities have been recorded in previous studies.
3. To be thorough, you check whether random draws of cysteines overlap with the feature you're curious about.
4. You calculate a p-value to quantify your confidence in the difference between these comparisons.
5. You repeat this procedure for every feature you're interested in
6. You Benjamini-Hochberg adjust your p-values, and ... et voilà!
