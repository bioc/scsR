\name{delColDf}
\alias{delColDf}
\title{
delColDf
}
\description{
Delete a specific column in the data frame.
}
\usage{
delColDf(df, colName)
}
\arguments{
    \item{df}{
      data frame
}
    \item{colName}{
      name of the column to be deleted (character vector)
}
}
\value{
input data frame with the column deleted.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	uuk_screen2 = delColDf(uuk_screen, "score")
}