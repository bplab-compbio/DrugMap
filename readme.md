# DrugMap

Welcome to the home-page of **DrugMap: A quantitative pan-cancer analysis of cysteine ligandability**! \n

<p align="center">
  <img src="https://github.com/bplab-compbio/DrugMap/blob/main/src/images/circos.png" width="700" height="600">
</p>

In order to systematically identify proteins in cancer which contain covalent opportunities for cysteine liganding, we used a mass-spectrometry based assay to profile the landscape of cysteine reactivity across 416 cancer models. 

<p align="center">
  <img src="https://github.com/bplab-compbio/DrugMap/blob/main/src/images/cysteine.architecture.png"  width="600" height="600">
</p>

Along the way, we began to reveal principles which determine cysteine liganding across diverse protein-structural and oncogenic contexts. More on this [here](https://www.cell.com/cell/abstract/S0092-8674(24)00318-0#secsectitle0020).

In order to make this resource maximally helpful to the wider community of cancer biologists, cysteine sleuths, and other interested travelers, we have released essential methodology underpinning our analyses, including:

1. The raw outputs of the engines we used to search our cysteine-targeted mass spectrometry data
2. The complete code we used to wrangle together the DrugMap data set at every level of normalization and analysis
3. The function and initial database for our tool Cysteine Set Enrichment Analysis (CSEA), as well as an example of its use
4. The code we used to train a neural network aimed at predicting cysteine ligandability
