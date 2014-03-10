\name{enrichment_geneSet}
\alias{enrichment_geneSet}
\title{
enrichment_geneSet
}
\description{
Computes the hypergeometric p-value that represents the enrichment of genesetA with genes of the genesetB.
}
\usage{
enrichment_geneSet(genesetA, genesetB, background=NULL, quiet=FALSE)
}
\arguments{
  \item{genesetA}{
    vector of sorted genes (vector of strings)
}
  \item{genesetB}{
    vector of sorted genes (vector of strings)
}
  \item{background}{
      vector of genes to be used as background (vector of strings)
}
  \item{quiet}{
    avoid print any message/warning (boolean)
}
}
\value{
the hypergeometric p-value that represents the enrichment of genesetA with genes of the genesetB. (integer)
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	data(uuk_screen_dh) 
	enrichment_geneSet(unique(uuk_screen$GeneID)[1:200], 
                   unique(uuk_screen_dh$GeneID)[1:400], 
                   intersect(uuk_screen$GeneID, uuk_screen_dh$GeneID))
}
