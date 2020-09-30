# Dincta

*Data INtegration and Cell Type Annotation of Single Cell Transcriptomes*

Check out the latest preprint of Dincta on [bioRxiv](https://www.biorxiv.org/content/early/2020/tobeappear)

# System requirements 

Dincta has been tested on R versions >= 3.5.0. Please consult the DESCRIPTION file for more details on required R packages. Dincta has been tested on OS X,  has not tested on Linux,  and Windows platforms.

# Installation

To run Dincta, open R and install directly from github using the following commands: 

```
library(devtools)
install_github("immunogenomics/dincta")
```

Installation may include compiling C++ code from source, so it can take a few minutes. 

# Usage/Demos

We made it easy to run Dincta in most common R analysis pipelines. 


## PCA matrix

The Dincta algorithm iteratively corrects PCA embeddings. To input your own low dimensional embeddings directly, set `do_pca=FALSE`. Dincta is packaged with a small dataset 

```
library(dincta)
my_dincta_embeddings <- DinctaMatrix(my_pca_embeddings, meta_data, "dataset", "cell_type", do_pca=FALSE)
```

## Normalized gene matrix

You can also run Dincta on a sparse matrix of library size normalized expression counts. Dincta will scale these counts, run PCA, and finally perform integration. 

```
library(dincta)
my_dincta_embeddings <- DinctaMatrix(normalized_counts, meta_data, "dataset", "cell_type")
```


## Dincta with two or more covariates

Dincta can integrate over multiple covariates. To do this, specify a vector covariates to integrate. 

```
my_dincta_embeddings <- DinctaMatrix(my_pca_embeddings, meta_data, c("dataset", "donor", "batch_id"), "cell_type", do_pca=FALSE)
```





