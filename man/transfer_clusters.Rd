\name{transfer_clusters}
\alias{transfer_clusters}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transfer seurat cluster labels from a SeuratObject
}
\description{
Given a SeuratObject, this function transfers cluster labels to a Visium tissue slide
}
\usage{
transfer_clusters(obj, "SeuratProject", plt1, "tissue_position_list.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seuratObj}{
%% SeuratObject
}
  \item{proj}{
%%  Project name. "SeuratProject" by default
}
  \item{plt}{
%%     DimPlot(obj, reduction="umap")
}
  \item{tissue_csv}{
%%    Standard space ranger pipeline output
}
}
