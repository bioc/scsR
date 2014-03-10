\name{launch_RSA}
\alias{launch_RSA}
\title{
launch_RSA
}
\description{
launch RSA sorting method
}
\usage{
launch_RSA(df, LB=-100, UB=100, reverse=FALSE, strScoreCol="", strGeneCol="Gene_ID", keepAllRSAReturnFields=FALSE)
}
\arguments{
  \item{df}{
    data frame containing the results of the siRNA experiment.
}
  \item{LB}{
    RSA lower bound (look KONIG paper). (number)
}
  \item{UB}{
    RSA upper bound (look KONIG paper). (number)
}
  \item{reverse}{
    whether to sort in ascending or descending order. (boolean)
}
\item{strScoreCol}{
    name of the column that contains the score of the screen (character vector)
}
\item{strGeneCol}{
    name of the column that contains the names of the genes in the screen (character vector)
}
\item{keepAllRSAReturnFields}{
    specify whether you want to keep all RSA columns in the output file. (boolean)
}
}
\value{
screen data frame with RSA sorting columns added.
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
	
	#extract the first 1000 lines in order to speed up the example
	screen = uuk_screen[1:1000,]
	
	screen_ranked <- launch_RSA(screen, strGeneCol="GeneID", strScoreCol="score")
}