\name{randomSortOnVal}
\alias{randomSortOnVal}
\title{
randomSortOnVal
}
\description{
randomize the order of the rows, on the values of a column
(e.g. randomized the rows, keeping close the rows having the same GeneID... i.e. sort the Genes of the screen in a random way).
}
\usage{
randomSortOnVal(screen, strColVal)
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{strColVal}{
    column of which the values have to be kept close to each other (character vector)
}
}
\value{
screen data frame sorted randomly on the defined column.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	
	# to speed up the example we use only the first 1000 rows
	uuk_screen_reduced = uuk_screen[1:1000,]

	screen = randomSortOnVal(uuk_screen_reduced, "GeneID")
}