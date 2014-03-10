\name{median_replicates}
\alias{median_replicates}
\title{
median_replicates
}
\description{
perform the median of the replicates (i.e. group by oligo sequence and takes the median of the score value ).
}
\usage{
median_replicates(screen, seedColName = "seed7", scoreColName = "score",
                                 geneColName = "GeneID", seqColName="siRNA_seq", spAvgColName = NULL)
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{seedColName}{
    name of the column that contains the seed of the sequence (character vector)
}
  \item{scoreColName}{
    name of the column that contains the score of the screen (character vector)
}
  \item{geneColName}{
    name of the column that contains the names of the genes in the screen (character vector)
}
\item{seqColName}{
    name of the column that contains the names of the sequences in the screen 
}
\item{spAvgColName}{
    name of the column that contains the names of the genes in the screen (character vector)
}
}
\value{
input data frame after having performed the median of the replicates
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	mr <- median_replicates(uuk_screen)
}
