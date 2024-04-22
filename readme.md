# DrugMap

Welcome to the home-page of **DrugMap: A quantitative pan-cancer analysis of cysteine ligandability**!  \\


<p align="center">
  <img src="https://github.com/bplab-compbio/DrugMap/blob/main/src/images/circos.png" width="700" height="600">
</p>

Toward systematically identifying proteins in cancer which contain covalent opportunities for cysteine liganding, we used a mass-spectrometry based assay to profile the landscape of cysteine reactivity across 416 cancer models. 

<p align="center">
  <img src="https://github.com/bplab-compbio/DrugMap/blob/main/src/images/cysteine.architecture.png"  width="600" height="600">
</p>

Along the way, we began to reveal principles which determine cysteine liganding across diverse protein-structural and oncogenic contexts. More on this [here](https://www.cell.com/cell/abstract/S0092-8674(24)00318-0#secsectitle0020).

To make this resource maximally helpful to the wider community of cancer biologists, cysteine sleuths, and other interested travelers, we have released essential methodology underpinning our analyses, including:

1. The [raw outputs](https://github.com/bplab-compbio/DrugMap/tree/main/Data) of the software that we used to search our cysteine-targeted mass spectrometry data
2. The complete, detailed [workflow](https://github.com/bplab-compbio/DrugMap/blob/main/src/misc/wrangle.m) that we used to wrangle together the DrugMap data set
3. The function and initial database for our tool [Cysteine Set Enrichment Analysis (CSEA)](https://github.com/bplab-compbio/DrugMap/tree/main/CSEA), as well as an example of its use
4. The [code](https://github.com/bplab-compbio/DrugMap/tree/main/src/neuralnet) we used to train a neural network aimed at predicting cysteine ligandability
