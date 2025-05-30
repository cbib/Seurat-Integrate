% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metrics_ari.R
\name{ScoreARI}
\alias{ScoreARI}
\alias{AddScoreARI}
\title{Score a clustering result with adjusted rand index}
\usage{
ScoreARI(object, cell.var, clust.var = "seurat_clusters")

AddScoreARI(object, integration, cell.var, clust.var = "seurat_clusters")
}
\arguments{
\item{object}{A Seurat object}

\item{cell.var}{The name(s) of the column(s) with cell type label variable
(must be in the object metadata). Multiple column names are accepted}

\item{clust.var}{The name of the column with cluster id assignment for each
cell (must be in the object metadata). Only one column name is accepted}

\item{integration}{name of the integration to score}
}
\value{
\code{ScoreARI}: a named array with as many values as there are
common strings between cell.var and the column names of the object's
metadata. Names are cell.var and values are ARI.

\code{AddScoreARI}: the updated Seurat \code{object} with the ARI score(s)
set for the integration.
}
\description{
Compute a score based on adjusted rand index between a clustering
result and one or more cell type label variable(s). 0 and 1 reflect a random
clustering and a perfect clustering as compared to cell type labelling
respectively.
}
\details{
ARI is rand index corrected for chance:
\deqn{\displaystyle ARI = \frac{RI - RI_{expected}}{max(RI) - RI_{expected}}}
More precisely, a contingency table is computed with the two variables
\eqn{L} and \eqn{C} of \eqn{r} and \eqn{s} elements respectively. For
\eqn{i \in [\![1,r]\!]} and \eqn{j \in [\![1,s]\!]}, \eqn{n_{ij}} is the
number of common samples (i.e. cells) between \eqn{L_i} and \eqn{C_j},
\eqn{a_i} is the number of samples in \eqn{L_i} and \eqn{b_j} is the number
of samples in \eqn{C_j}. The ARI is:
\deqn{\displaystyle ARI = \frac{\left. \sum_{ij} \binom{n_{ij}}{2} - \left(\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2}\right) \right/ \binom{n}{2} }{ \left. \frac{1}{2} \left(\sum_i \binom{a_i}{2} + \sum_j \binom{b_j}{2}\right) - \left(\sum_i \binom{a_i}{2} \sum_j \binom{b_j}{2}\right) \right/ \binom{n}{2}}}
}
\note{
The metric is symmetric. Switching cell.var with clust.var will return
the same value.
}
\references{
Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A.,
Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M.,
Colomé-Tatché, M. & Theis, F. J. Benchmarking atlas-level data integration in
single-cell genomics. Nat Methods 19, 41–50 (2021).
\href{https://doi.org/10.1038/s41592-021-01336-8}{DOI}
}
