\name{read_h5}
\alias{read_h5}
\title{Takes h5 file and creates a Seurat object, normalizes, does variable feature selection, scales the data, and performs PCA.}
\usage{
get_expression(seuratObj, proj, feature, tissue_csv)
}
\arguments{
\item{h5}{Seurat object.}

\item{proj}{Project name}

\item{a}{Assay}

}
\value{
Seurat object
}
\description{
Takes h5 file and creates a Seurat object, normalizes, does variable feature selection, scales the data, and performs PCA.
}
\examples{
read_h5("example.h5","my_project","RNA")
}
