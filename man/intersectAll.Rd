\name{intersectAll}
\alias{intersectAll}
\title{
intersectAll
}
\description{
  intersect several vectors that can be passed as arguments of the functino
}
\usage{
intersectAll(...)
}
\arguments{
  \item{...}{
    vectors to intersect
}
}
\value{
vector that results from the intersection of the input vectors
}
\author{
Andrea Franceschini
}
\examples{
	intersectAll(c(1,2,3,4), c(1,2), c(2,3,4))
}
