\name{bydf}
\alias{bydf}
\alias{bydfa}
\title{
bydf
}
\description{
apply a function to a group of rows in the input data frame (similar to the sql group by statememnt).
}

\details{
  The methods currently depend on the type of \code{to}:
  \describe{
	\item{bydf}{
		apply a function to a group of rows in the input data frame (similar to the sql group by statememnt).
Put the results of this function in a new data frame that is returned as output.
	}
	\item{bydfa}{
		apply a function to a group of rows in the input data frame (similar to the sql group by statememnt).
return the same data frame with an additional column with the results of the function.	}
}
}

\usage{

bydf(df, groupColName, valColName, fun, newColName="temp_by_col_name")

bydfa(df, groupColName, valColName, fun, newColName="temp_by_col_name")



}
\arguments{
  \item{df}{
input data frame
}
  \item{groupColName}{
    name of the column to be used for grouping the rows (character vector)
}
  \item{valColName}{
    name of the column containing the values to be inserted in the function (character vector)
}
  \item{fun}{
    function to be applied (function)
}
\item{newColName}{
    name of the colum that contains the result of the function (character vector)
}
}
\value{
bydf: data frame with the function applied to the grouping 
bydfa: input data frame with an additional column with the results of the function applied to the grouping.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	screen=add_seed(uuk_screen[1:1000,])
	screen_sd = bydf(screen, groupColName="seed7", "score", sd, "sd")
}
