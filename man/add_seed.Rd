\name{add_seed}
\alias{add_seed}
\title{
add_seed
}
\description{
This method takes in input a dataframe containing the results of an siRNA screen. This screen must contain the siRNA sequences in a dedicated column (the sequences have to be provided in the guide/antisense orientation).
Then it adds a column with the seed of the siRNA sequences.
}
\usage{
add_seed(df, seqColName="siRNA_seq", seedLength=7, startPosition=2)
}
\arguments{
  \item{df}{
    Dataframe containing the results of the siRNA screen.
}
  \item{seqColName}{
    character vector with the name of the column that contains the siRNA sequences (the sequences have to be provided in the guide/antisense orientation).
}
  \item{seedLength}{
    length of the seed in nucleotides (by default 7 bases)  (integer)
}
  \item{startPosition}{
    position in the siRNA sequence where the seed starts (by default position 2) (integer))
}
}
\value{
screen data frame with the seed column added.
}
\author{
Andrea Franceschini
}
\examples{
 data(uuk_screen)
 seed_uuk_screen = add_seed(uuk_screen[1:100,])
}