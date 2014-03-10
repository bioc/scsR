\name{replace_non_null_elements}
\alias{replace_non_null_elements}
\title{
replace_non_null_elements
}
\description{
replace the element of the input vector with the element of the replacementVector (whenever these elements are not empy/null)
}
\usage{
replace_non_null_elements(inputVect, replacementVect)
}
\arguments{
  \item{inputVect}{
    data frame containing the results of the siRNA experiment.
}
  \item{replacementVect}{
	replacement vector (vector)
}
}
\value{
input vector with the replaced values.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 1000 rows
	uuk_screen_reduced = uuk_screen[1:1000,]

	# replace all the scores with 1, except the first 100 scores of the vector
	nv <- replace_non_null_elements(uuk_screen_reduced$score, c(rep(NA, 100), rep(1, nrow(uuk_screen_reduced)-100) ))
}
