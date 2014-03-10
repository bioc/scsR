\name{create_sd_matrix}
\alias{create_sd_matrix}
\title{
create_sd_matrix
}
\description{
We observed that the standard deviation of the oligos that share the same seed do change relative to their average score.
In principle we could plot this information on a graph (x-axis = average of the oligos that share the same seed, y-axes = standard deviation of the oligos).
We do provide this utility method to condense this information in a matrix (that reports the quantiles of the standard deviation for every score interval).
}
\usage{
create_sd_matrix(screen, seedColName="seed7", scoreColName="score")
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{seedColName}{
    name of the column that contains the seeds sequences of the screen (character vector)
}
  \item{scoreColName}{
    name of the column that contains the score of the screen (character vector)
}
}
\value{
matrix that reports the quantiles of the standard deviation for every score interval.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	
	
	# to speed up the example we use only the first 100 rows
	uuk_screen_reduced = uuk_screen[1:100,]

	screen = add_seed(uuk_screen_reduced)
	sd_matrix = create_sd_matrix(screen)
}