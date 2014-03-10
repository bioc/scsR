\name{check_consistency}
\alias{check_consistency}
\title{
check_consistency
}
\description{
This method takes an siRNA screen as input and check its consistency (i.e. check that the format of the data is suitable for the usage with our scsR package).
The method prints meaningful warnings for every inconsistency that can be detected
}
\usage{

check_consistency(screen, scoreColName = "score", geneColName = "GeneID",
                                  seqColName="siRNA_seq")
}
\arguments{
  \item{screen}{
Dataframe containing the results of the siRNA scree
}
  \item{scoreColName}{
    name of the column that contains the score of the screen (character vector)
}
  \item{geneColName}{
    name of the column that contains the gene identifier of the screen (character vector)
}
  \item{seqColName}{
    name of the column that contains the siRNA sequences in the screen. (the sequences have to be provided in the guide/antisense orientation and each sequence must be in the format of a character vector, i.e. a simple string).  (character vector)
}
}
\value{
return the data frame passed in input with possible consistency corrections.
}
\author{
Andrea Franceschini
}
\examples{
data(uuk_screen)
uuk_screen <- check_consistency(uuk_screen)
}
