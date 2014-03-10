\name{add_rank_col}
\alias{add_rank_col}
\title{
add_rank_col
}
\description{
This method takes in input a dataframe containing the results of an siRNA screen.
Then it adds a set of column that are useful for sorting to the dataframe.
At the moment the following sorting columns are provided: 
- column with the median value of the siRNA score for each gene
- columns that comes out from the execution of the RSA sorting method (Renate Konig et al.)
}
\usage{
add_rank_col(screen, reverse=FALSE,  scoreColName="score", geneColName="GeneID")
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{reverse}{
    boolean specifying the direction of the sorting (from the lowest scores to the highest score or vice versa)
}
  \item{scoreColName}{
    character vector with the name of the column that contains the score of the screen
}
  \item{geneColName}{
    character vector withname of the column that contains the names of the genes in the screen
}
}
\value{
screen data frame with sorting columns added.
}
\references{
A probability-based approach for the analysis of large-scale RNAi screens. 
Renate Konig et al. Nature Methods 2007
}
\author{
Andrea Franceschini
}
\examples{
 data(uuk_screen)
 uuk_screen_ranked = add_rank_col(uuk_screen[1:100,])
}