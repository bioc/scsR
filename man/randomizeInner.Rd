\name{randomizeInner}
\alias{randomizeInner}
\title{
randomizeInner
}
\description{
randomize an inner field (e.g. the scores of the oligos of a gene), keeping unaltered the order of the outer field (e.g. the genes)
}
\usage{
randomizeInner(df, baseColStr, sortColStr, reverse = FALSE) 
}
\arguments{
    \item{df}{
    input data frame
}
\item{baseColStr}{
    name of the column that represents the outer field (e.g. the genes) (character vector)
}
\item{sortColStr}{
    name of the column that represents the inner field (e.g. the scores of the oligos of a gene) (character vector)
}
\item{reverse}{
    specify the direction of the sorting (boolean)
}
}
\value{
data frame with the randomized rows.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 1000 rows
	uuk_screen_reduced = uuk_screen[1:1000,]

	screen <- randomizeInner(arrange(uuk_screen_reduced, GeneID), "GeneID", "score")
}
