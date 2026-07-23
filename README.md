# SeuratIntegrate <!-- omit in toc -->
R package expanding integrative analysis capabilities of Seurat by providing seamless access to popular integration methods. It also implements an integration benchmarking toolkit that gathers well-established performance metrics to help select the most appropriate integration.

Examples, documentation, memos, etc. are available on the SeuratIntegrate's [website](https://cbib.github.io/Seurat-Integrate/).

SeuratIntegrate provides support to R- and Python-based integration methods. The table below summarizes which methods are compatible with SeuratIntegrate:

<table id="table1" style="margin: 0px auto;"><thead>
  <caption>Table 1: Supported integration methods</caption>
  <tr>
    <th></th>
    <th>Package</th>
    <th>Method</th>
    <th>Function</th>
  </tr></thead>
<tbody>
  <tr>
    <td rowspan="6">R</td>
    <td rowspan="3">SeuratIntegrate</td>
    <td><a href="https://www.bioconductor.org/packages/release/bioc/html/sva.html" target="_blank" rel="noopener noreferrer">ComBat</a></td>
    <td><code>CombatIntegration()</code></td>
  </tr>
  <tr>
    <td><a href="https://cran.r-project.org/web/packages/harmony/index.html" target="_blank" rel="noopener noreferrer">Harmony</a></td>
    <td><code>HarmonyIntegration()</code></td>
  </tr>
  <tr>
    <td><a href="https://www.bioconductor.org/packages/release/bioc/html/batchelor.html" target="_blank" rel="noopener noreferrer">MNN</a></td>
    <td><code>MNNIntegration()</code></td>
  </tr>
  <tr>
    <td rowspan="2">Seurat</td>
    <td><a href="https://cran.r-project.org/web/packages/Seurat/index.html" target="_blank" rel="noopener noreferrer">CCA</a></td>
    <td><code>CCAIntegration()</code></td>
  </tr>
  <tr>
    <td><a href="https://cran.r-project.org/web/packages/Seurat/index.html" target="_blank" rel="noopener noreferrer">RPCA</a></td>
    <td><code>RPCAIntegration()</code></td>
  </tr>
  <tr>
    <td>SeuratWrappers</td>
    <td><a href="https://github.com/satijalab/seurat-wrappers" target="_blank" rel="noopener noreferrer">FastMNN</a><br>(<a href="https://bioconductor.org/packages/release/bioc/html/batchelor.html" target="_blank" rel="noopener noreferrer">batchelor</a>)</td>
    <td><code>FastMNNIntegration()</code></td>
  </tr>
  <tr>
    <td rowspan="5">Python</td>
    <td rowspan="5">SeuratIntegrate</td>
    <td><a href="https://github.com/Teichlab/bbknn" target="_blank" rel="noopener noreferrer">BBKNN</a></td>
    <td><code>bbknnIntegration()</code></td>
  </tr>
  <tr>
    <td><a href="https://github.com/scverse/scvi-tools" target="_blank" rel="noopener noreferrer">scVI</a></td>
    <td><code>scVIIntegration()</code></td>
  </tr>
  <tr>
    <td><a href="https://github.com/scverse/scvi-tools" target="_blank" rel="noopener noreferrer">scANVI</a></td>
    <td><code>scANVIIntegration()</code></td>
  </tr>
  <tr>
    <td><a href="https://github.com/brianhie/scanorama" target="_blank" rel="noopener noreferrer">Scanorama</a></td>
    <td><code>ScanoramaIntegration()</code></td>
  </tr>
  <tr>
    <td><a href="https://github.com/theislab/scarches" target="_blank" rel="noopener noreferrer">trVAE</a></td>
    <td><code>trVAEIntegration()</code></td>
  </tr>
</tbody></table>


- [Installation](#installation)
- [Preparations](#preparations)
  - [Setup Python environments](#setup-python-environments)
  - [Setup a `SeuratObject`](#setup-a-seuratobject)
  - [Facultative dependencies](#facultative-dependencies)
- [SeuratIntegrate usage](#seuratintegrate-usage)
  - [Integrate datasets](#integrate-datasets)
  - [Post-process integration outputs](#post-process-integration-outputs)
  - [Compare integrations](#compare-integrations)
- [Getting help and advice](#getting-help-and-advice)
- [Citing](#citing)


## Installation

First, use the yaml files to creates environments using conda
https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html


Install SeuratIntegrate from github directly:
```R
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!require("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("cbib/Seurat-Integrate", dependencies = NA, repos = BiocManager::repositories()) 
```

## Preparations

### Setup Python environments

To use Python methods, run the following commands (once) to set up the necessary conda environments:

```R
library(SeuratIntegrate)

# Create envrionments
UpdateEnvCache("bbknn")
UpdateEnvCache("scvi")     # also scANVI
UpdateEnvCache("scanorama")
UpdateEnvCache("trvae")

# Show cached environments
getCache()
```

Environments are persistently stored in the cache and the **`UpdateEnvCache()` commands should not need to be executed again**.

While these environments should work well in most cases, conda's dependencies occasionally encounter conflicts. Manual adjustment might be needed. You may find helpful information in [this vignette](https://cbib.github.io/Seurat-Integrate/articles/setup_and_tips.html#troubleshouting-with-conda).

### Setup a `SeuratObject`

To integrate data with SeuratIntegrate, you need to preprocess your `SeuratObject` until you obtain at least a PCA. Importantly, the `SeuratObject` must have its **layers split by batches**.

<details>
  <summary>Not familiar with Seurat?</summary>
  
 Have a look at Seurat's [website](https://satijalab.org/seurat/articles/get_started_v5_new), especially the tutorials covering [SCTransform](https://satijalab.org/seurat/articles/sctransform_vignette) and [integrative analyses](https://satijalab.org/seurat/articles/seurat5_integration).
</details>
<br>

To fully benefit from the benchmarking toolkit, you'll need **cell-type annotations of sufficient quality** to be considered suitable as ground truth.

### Facultative dependencies

The benchmarking toolkit can benefit from additional dependencies:
```R
# required to test for k-nearest neighbour batch effects
remotes::install_github('theislab/kBET')

# fast distance computation
install.packages('distances')

# faster Local Inverse Simpson’s Index computation
remotes::install_github('immunogenomics/lisi')
```

## SeuratIntegrate usage
### Integrate datasets

When your `SeuratObject` is ready, you can launch multiple integrations (from [Table 1](#table1)) with a single command. `DoIntegrate()` provides a flexible interface to customise integration-specific parameters and to control over associated data and features.

```R
seu <- DoIntegrate(seu,
       # ... integrations
         CombatIntegration(layers = "data"),
         HarmonyIntegration(orig = "pca", dims = 1:30),
         ScanoramaIntegration(ncores = 4L, layers = "data"),
         scVIIntegration(layers = "counts", features = Features(seu)),
       # ...
       use.hvg = TRUE,    # `VariableFeatures()`
       use.future = c(FALSE, FALSE, TRUE, TRUE)
)
```

In this example, all integration methods will use the variable features as input, with the exception of `scVIIntegration()` which is set to use all features (`features = Features(seu)`). `CombatIntegration()` will correct the normalised counts (`layers = "data"`), while `scVIIntegration()` will train on raw counts (`layers = "counts"`).

`use.future` must be **`TRUE` for Python methods**, and `FALSE` for R methods (see [Table 1](#table1)). 

### Post-process integration outputs

Integration methods produce one or several outputs. Because they can be of different types, the following table indicates the post-processing steps to generate a UMAP.

<a id="table2"></a>
<caption>Table 2: Output types and processing</caption>

| **Output type**       | **Object name** |                             **Processing** |
|-----------------------|:---------------:|-------------------------------------------:|
| Corrected counts      |     `Assay`     |   `ScaleData()` ➔ `RunPCA()` ➔ `RunUMAP()` |
| Dimensional reduction |    `DimReduc`   |                                `RunUMAP()` |
| KNN graph             |     `Graph`     |      `RunUMAP(umap.method = "umap-learn")` |

Output types are summarized for each method in the [Memo vignette](https://cbib.github.io/Seurat-Integrate/articles/memo_integration.html) about integration methods


### Compare integrations

SeuratIntegrate incorporates 11 scoring metrics: 6 quantify the degree of batch mixing (*batch correction*), while 5 assess the preservation of biological differences (*bio-conservation*) based on ground truth cell type labels.

To score your integrations, you must process their outputs as in the **Processing** column of [Table 2](#table2). You'll also need to get a graph by running `FindNeighbors(return.neighbor = TRUE)` (this [vignette](https://cbib.github.io/Seurat-Integrate/articles/SeuratIntegrate.html#post-processing) provides further guidance).

Then, scores can be obtained using the function `Score[score_name]()`, or directly saved in the Seurat object using the `AddScore[score_name]()` as follows:

```R
# save the score in a variable
rpca_score <- ScoreRegressPC(seu, reduction = "[dimension_reduction]")  #e.g. "pca"

# or save the score in the Seurat object
seu <- AddScoreRegressPC(seu, integration = "[name_of_integration]", reduction = "[dimension_reduction]")
```
It is worth noting that the **unintegrated version must also be scored** to perform a complete comparative analysis. When scores have been computed, they can be used to compare the integration outputs. See this [vignette](https://cbib.github.io/Seurat-Integrate/articles/memo_score.html) for a complete overview of available scores.

The advantage of the `AddScore` over the `Score` functions is that they facilitate score scaling and plotting:

```R
# scale
seu <- ScaleScores(seu)

# plot
PlotScores(seu)
```

## Getting help and advice

Examples, documentation, memos, etc. are available on SeuratIntegrate's [website](https://cbib.github.io/Seurat-Integrate/).

If you encounter a bug, please create an [issue on GitHub](https://github.com/cbib/Seurat-Integrate/issues). Likewise if you have a specific comment or question not covered on the website.

## Citing

If you find SeuratIntegrate useful, please consider citing:

> Specque, F., Barré, A., Nikolski, M., & Chalopin, D. (2025). SeuratIntegrate: an R package to facilitate the use of integration methods with Seurat. *Bioinformatics*. doi: [10.1093/bioinformatics/btaf358](https://doi.org/10.1093/bioinformatics/btaf358)
