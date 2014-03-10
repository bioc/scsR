\name{split_df}
\alias{split_df}
\title{
split_df
}
\description{
You can use this function to extract only certain rows from the data frame.
Give a value in the strIdCol, we include in the output data frame only the rows specified in the linesToGet vector
(e.g. for every gene extract only the first two oligos)
}
\usage{
split_df(df, strIdCol, linesToGet)
}
\arguments{
 \item{df}{
    input data frame
}
\item{strIdCol}{
    name of the column containing the identifiers of the groups (character vector)
}
\item{linesToGet}{
    vector containing the numbers of the lines to retrieve (vector)
}

}
\value{
subset of the input data frame (only the rows requested are included in the subset)
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 1000 rows
	uuk_screen_reduced = uuk_screen[1:1000,]


	uuk_screen_firstOligo <- split_df(uuk_screen_reduced, "GeneID", c(1))
}
