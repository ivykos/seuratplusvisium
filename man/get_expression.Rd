\name{get_expression}
\alias{get_expression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize gene expression spatially
}
\description{
This function visualizes the expression of a gene of interest on a Visium tissue slide. Note this does not yet
work for SeuratObjects made of merged datasets.
}
\usage{
get_expression(obj, "SeuratProject", "Cx3cr1", "tissue_position_list.csv")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{seuratObj}{
%%     SeuratObject
}
  \item{proj}{
%%     Project name (SeuratProject by default)
}
  \item{feature}{
%%     Gene or feature of interest
}
  \item{tissue_csv}{
%%    Standard space ranger output
}
}
