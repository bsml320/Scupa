# Scupa

Scupa is an R package for immune cell polarization analysis of scRNA-seq data.

Scupa relies on Universial Cell Embedding (UCE) to measure the polarization of 14 immune cell types:

B, NK, CD8T, CD4T, Treg, Tgd, pDC, cDC1, cDC2, MigDC, LC (Langerhans cell), Macrophage, Monocyte, Neutrophil.

With UCE's property, Scupa can be applied to various animals for unified polarization analysis.

![Scupa design and immune cell polarization states](inst/figure/scupa.png)

## Installation

```r
devtools::install_github("bsml320/scupa")
```

## Usage and vignettes

### Prepare UCE (using python)

See the vignette [vignette_run_UCE_ifnb.ipynb](inst/notebook/vignette_run_UCE_ifnb.ipynb) for a quick tutorial.

Before running Scupa, it is necessary to generate the UCE embeddings for the input scRNA-seq dataset. Please install UCE (https://github.com/snap-stanford/UCE) and run it on the input dataset.

### Run Scupa (using R)

See the vignette [vignette_scupa_ifnb.ipynb](inst/notebook/vignette_scupa_ifnb.ipynb) for a quick tutorial.

The R package schard (https://github.com/cellgeni/schard/) helps convert the UCE output h5ad data to a Seurat object. Otherwise, the users could create a Seurat object with UCE embeddings saved in an assay or dimensional reduction structure on their own.

```r
seuobj <- schard::h5ad2seurat('output_uce_adata.h5ad')
```

After conversion, run *MeasurePolar* function to analyze the polarization of input cells to different polarization states.

```r
library(scupa)
# Replace input_type with one of: 
#   B, NK, CD8T, CD4T, Treg, Tgd, pDC, cDC1, cDC2, MigDC, LC, Macro, Mono, Neu.
# Replace cell_names with the names of unpolarized cells, or NA if uncertain.
#   See ?MeasurePolar for details
seuobj <- MeasurePolar(seuobj, celltype=input_type, embedding='uce', unpolarized_cell=cell_names)
```

## Output

By default, Scupa outputs a Seurat object with the updated metadata. Each cell type has 4-6 polarization states. For each polarization state, the polarization scores (range 0-1), p-values, and adjusted p-values are calculated for all cells. Larger scores indicate stronger polarization.

The scores and p-values can be visualized by Seurat functions FeaturePlot, VlnPlot, etc.

## Citations

If you use Scupa, please consider citing following papers:

* Liu W, Zhao Z. Scupa: immune cell polarization analysis using the single-cell foundation model. BioRxiv. 2024. doi: comming soon.
* Cui A, Huang T, Li S, et al. Dictionary of immune responses to cytokines at single-cell resolution. Nature. 2024;625(7994):377-384. doi:10.1038/s41586-023-06816-9
* Rosen Y, Roohani Y, et al. Universal Cell Embeddings: A Foundation Model for Cell Biology. BioRxiv. 2023. doi: 10.1101/2023.11.28.568918

