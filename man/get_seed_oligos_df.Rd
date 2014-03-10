\name{get_seed_oligos_df}
\alias{get_seed_oligos_df}
\title{
get_seed_oligos_df
}
\description{
  This function returns the screen, that is given in input, with additional columns about the possible off-targets/seed effect of each oligos.
  The seed effect is computed excluding the current oligo.
}
\usage{
get_seed_oligos_df(screen, seedColName="seed7",  scoreColName="score", geneColName="GeneID", gene_interval = c(1,100), 
                                 min_oligos_x_gene=4, min_oligos_x_statistics=4, random=FALSE, kolmogorovSampleSize=5000, progress_bar=FALSE )
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment (sorted by significance).
}
  \item{seedColName}{
    specify the direction of the sorting (from the lowest scores to the highest score or vice versa) (character vector)
}
  \item{scoreColName}{
    name of the column that contains the score of the screen	 (character vector)
}
  \item{geneColName}{
    name of the column that contains the names of the genes in the screen	 (character vector)
}
  \item{gene_interval}{
    apply the analysis only to the genes that are included in this interval (the screen must be sorted by significance and the interval has to be intended from the best hits to the worst hits).	(vector of integer)
}
  \item{min_oligos_x_gene}{
    minimum number of oligos that a gene must have in order to be included in the analysis (integer)
}
  \item{min_oligos_x_statistics}{
    minimum number of oligos with the same seed that is required in order to apply a statistics (otherwise 0 is returned). (integer)
}
  \item{random}{
    randomize the genes of the screen (boolean)
}
  \item{progress_bar}{
    print progress bar (boolean)
}
  \item{kolmogorovSampleSize}{
    sample size to be used for the Kolmogorov Smirnov statistics (i.e. the number of genes that we consider to be enough in order to infer the correct distribution of the genome-wide screen. The higher this number, the slower the computation).
    If this variable is left to NULL the Kolmogorov statistics is disabled (integer)
}
}
\value{
screen, that is given in input, with additional columns about the possible off-targets/seed effect of each oligos. (data frame)
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 100 rows
	uuk_screen_reduced = uuk_screen[1:1000,]

	uuk_screen <- add_seed(uuk_screen_reduced)
	sodf = get_seed_oligos_df(uuk_screen)
}

