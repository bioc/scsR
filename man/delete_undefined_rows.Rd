\name{delete_undefined_rows}
\alias{delete_undefined_rows}
\title{
delete_undefined_rows
}
\description{
method to delete the rows that contain undefined values in some specific columns.
}
\usage{
delete_undefined_rows(df, colNames, quiet=FALSE)
}
\arguments{
  \item{df}{
    data frame
}
  \item{colNames}{
    vector with the names of the column that must be defined (i.e. their values cannot be NULL, NA, NaN or zero-length strings) (vector of strings)
}
  \item{quiet}{
    specify whether to avoid printing warnings. (boolean)
}
}
\value{
data frame without the rows that contain at least one undefined value in the column list
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	screen <- delete_undefined_rows(uuk_screen, colNames=c("score", "GeneID"))
}