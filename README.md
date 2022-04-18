# seuratplusvisium

Seuratplusvisium is a set of R functions built on top of Seurat and designed to visualize cell cluster labels from Seurat's guided clustering workflow on a Visium
tissue slide.

Alongside this functionality, it can also visualize expression of any gene spatially.

The easiest way to install seuratplusvisium is to use the ```install_github()``` function from the devtools R package. Simply copy and paste the following code snippet to your R session.

```
library(devtools)
devtools::install_github('ivykos/seuratplusvisium')
library(seuratplusvisium)
```
Note: This package is still in early stages of development. It will be refined and more features will be added over time. 
