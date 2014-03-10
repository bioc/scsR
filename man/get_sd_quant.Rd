\name{get_sd_quant}
\alias{get_sd_quant}
\title{
get_sd_quant
}
\description{
This method scan the quantile standard deviation matrix (produced by create-sd-matrix function) and finds the quantile of the given standard deviation and average score
}
\usage{
get_sd_quant(sdval, score, sd_matrix)
}
\arguments{
  \item{sdval}{
    standard deviation (number)
}
  \item{score}{
    average score (number)
}
  \item{sd_matrix}{
    standard deviation quantile matrix (matrix)
}
}
\value{
number from 1 to 20 that represents the quantile of the standard deviation in the given score range (1 corresponds to 0.05 percent). (integer)
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 2500 rows
	uuk_screen_reduced = uuk_screen[1:2500,]

	screen = add_seed(uuk_screen_reduced)
	sd_matrix = create_sd_matrix(screen)
	quant <- get_sd_quant(0.3, 0.9, sd_matrix)
}