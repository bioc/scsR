\name{renameColDf}
\alias{renameColDf}
\title{
renameColDf
}
\description{
rename the column of a data frame
}
\usage{
renameColDf(df, colOldName, colNewName)
}
\arguments{
  \item{df}{
    input data frame
}
  \item{colOldName}{
    name of the column that has to be changed (character vector)
}
  \item{colNewName}{
    new column name (character vector)
}
}
\value{
input data frame with the name of the column changed
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	screen <- renameColDf(uuk_screen, "score", "my_z_score")
}
