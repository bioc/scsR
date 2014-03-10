\name{compare_sorted_geneSets}
\alias{compare_sorted_geneSets}
\title{
compare_sorted_geneSets
}
\description{
This method can be used to compare the performances of two different sorted gene vectors (A1 and A2) relative to a reference vector (B).
To perform the comparison we use n best hits from genesetA1 and genesetA2. n is defined as the number of elements of the smallest of the two vectors(after intersecting it with the background).
For the comparison see also the enrichment_geneSet method.
}
\usage{
compare_sorted_geneSets(genesetA1, genesetA2, genesetB, background, limA=NULL, limB=NULL)
}
\arguments{
  \item{genesetA1}{
    vector of sorted genes (character vector)
}
  \item{genesetA2}{
    vector of sorted genes (character vector)
}
  \item{genesetB}{
    vector of genes to be used as reference (character vector)
}
  \item{background}{
    vector of genes to be used as background (character vector)
}
  \item{limA}{
    limit the number of genes of the vector genesetA1 to the first limA genes (integer)
}
  \item{limB}{
    limit the number of genes of the vector genesetB1 to the first limB genes (integer)
}
}
\author{
Andrea Franceschini
}
\examples{
data(uuk_screen)
data(uuk_screen_dh) 
compare_sorted_geneSets(unique(uuk_screen$GeneID)[1:200], 
                        unique(arrange(add_rank_col(uuk_screen), log_pval_rsa)$GeneID)[1:200], 
                              unique(uuk_screen_dh$GeneID)[1:400], 
                              intersect(uuk_screen$GeneID, uuk_screen_dh$GeneID)
                              )
}
