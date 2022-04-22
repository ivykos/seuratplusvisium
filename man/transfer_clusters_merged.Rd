\name{transfer_clusters_merged}
\alias{transfer_clusters_merged}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transfer seurat cluster labels from a merged dataset
}
\description{
Given a SeuratObject consisting of merged datasets, this function transfers cluster labels to a Visium tissue slide
}
\usage{
transfer_clusters_merged(obj, "CTRL", "tissue_position_list.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seuratObj}{
%% SeuratObject
}
  \item{data_set}{
%%     metadata added to SeuratObject to indicate which data set you're interested in
}
  \item{tissue_csv}{
%%     Standard output from Space Ranger pipeline
}
}

