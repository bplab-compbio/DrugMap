Here, we used machine learning to probe a fundamental question in drug discovery... **can structural information be used to infer cysteine ligandability**? Where do we start?

<p align="center">
  <img src="https://github.com/bplab-compbio/DrugMap/blob/main/src/images/structural.mapping.png" width="700" height="300">
</p>

We first aligned as many cysteines to the [Protein Data Bank](https://www.rcsb.org/) as possible and deconvolved protein structures into their corresponding rotamers, when available. For each cysteine, we chose representative structures while optimizing for both 1.) completeness of structural coverage and 2.) resolution.

<p align="center">
  <img src="https://github.com/bplab-compbio/DrugMap/blob/main/src/images/neural.net.png" width="400" height="300">
</p>

We then fed a deep network both geometric and vectorized data, all encapsulating unique dimensions of a cysteine's structural locale, to allow the network to learn whether a cysteine is ligandable or not. We hope that users will find the structural data deposited herein to be a good starting point to train more sophisticated models in the future!
