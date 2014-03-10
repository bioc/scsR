\name{plot_screen_hits}
\alias{plot_screen_hits}
\title{
plot_screen_hits
}
\description{
Gene-seed plot:
plot the genes of the siRNA screen (x-axis) together with a representation of the effect of the seed of their oligos (circles).
The position on the y-axis of the circles refers to the average score of the oligos of the gene that share the same seed.
The dimension of the circles refers to the number of oligos that share the same seed in the screen (the higher the number of oligos with the same seed, the bigger is the circle).
This graph can be used to look by eyes at the effect of the seeds on the genes.
}
\usage{
plot_screen_hits(screen, output_file=NULL, geneScoreColName="median", seedColName="seed7",
                      scoreColName="score", geneColName="GeneID", gene_interval = c(1,100), 
                      min_oligos_x_gene=4, min_oligos_x_statistics=4, random=FALSE, kolmogorovSampleSize=5000, 
                      ylab="score", xlab="gene", ylim=c(-4,4), graph_highest_count_thr=16, progress_bar=FALSE)
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{output_file}{
    specify the direction of the sorting (from the lowest scores to the highest score or vice versa) (character vector)
}
  \item{scoreColName}{
    name of the column that contains the score of the screen (character vector)
}
  \item{geneColName}{
    name of the column that contains the names of the genes in the screen (character vector)
}
  \item{seedColName}{
    name of the column that contains the seeds of the siRNA sequences in the screen (character vector). 
(the sequences have to be provided in the guide/antisense orientation and each sequence must be in the format of a character vector, i.e. a simple string)
}
  \item{geneScoreColName}{
    name of the column that contains the global scores of the genes (i.e. the column that contains the median or the average value of the oligos ) (character vector)
}
  \item{gene_interval}{
    display in the graph only the genes that are included in this interval (the screen must be sorted by significance and the interval has to be intended from the best hits to the worst hits). (vector)
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
  \item{kolmogorovSampleSize}{
    sample size to be used for the Kolmogorov Smirnov statistics (i.e. the number of genes that we consider to be enough in order to infer the correct distribution of the genome-wide screen. The higher this number, the slower the computation).
    If this variable is left to NULL the Kolmogorov statistics is disabled (integer)
}

  \item{ylab}{
    label of the graph y-axis (character vector)
}
  \item{xlab}{
    label of the graph x-axis (character vector)
}
  \item{ylim}{
    ylim of the graph (vector)
}
  \item{graph_highest_count_thr}{
    maximum number of oligos to be used in order to display the largest circle in the graph (number)
}
  \item{progress_bar}{
    whether to show a progress bar or not (Boolean)
}
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 100 rows
	uuk_screen_reduced = uuk_screen[1:1000,]


	screen = add_rank_col(add_seed(uuk_screen_reduced))
	plot_screen_hits(screen)

	# The screen has to be sorted. In our case it is already sorted via median.
	# In order to sort the screen you can use our add_rank_col method
	# example: arrange(add_rank_col(screen), median) 
}

