\name{seed_removal}
\alias{seed_removal}
\title{
seed_removal
}
\description{
In certain cases we may want to select only the siRNA oligos that we are sure to NOT have strong detectable seed effect.
This is different than using our seed_correction method because we don't correct the scores based on an additivity assumption, but simply remove the oligos that shows a detectable off-target effect.
In principle such function should lead to the identification of few very reliable hits, but we will loose several potential hits (that can probably be detected using our seed_correction method). Hence we suggest the user to use first this seed_removal function, and then also our seed_correction method.
}
\usage{
seed_removal(screen,seedColName="seed7", scoreColName="score", geneColName="GeneID",
                            min_siRNAs_x_seed=4, remove_unrepresented_seeds=TRUE, lower_bound_threshold = -0.5,
                            higher_bound_threshold = 0.5, min_oligos_x_gene_threshold = 2, useMedian=FALSE, removeGenes=FALSE, include_current_gene=FALSE, progress_bar=FALSE)
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{seedColName}{
    name of the column that contains the seed of the siRNA oligo sequences. (character vector)
(the sequences have to be provided in the guide/antisense orientation and each sequence must be in the format of a character vector, i.e. a simple string)
}
  \item{scoreColName}{
    name of the column that contains the score of the screen (character vector)
}
  \item{geneColName}{
    name of the column that contains the names of the genes in the screen (character vector)
}
  \item{min_siRNAs_x_seed}{
    minimum number of oligos x seed that are required in order to execute the analysis (i.e. compute their average score and remove them if the score is outside the boundaries defined in the lower_bound_threshold and higher_bound_threshold parameters) (integer)
}
  \item{remove_unrepresented_seeds}{
    if set to TRUE, we remove all the seeds that are found in the screen in less than the number of oligos specified in the min_siRNAs_x_seed parameter. You should use this option if you think that it is not advisable to rely on siRNAs having seeds that are present in only few oligos, because we are not able to estimate precisely their seed effect and hence we cannot detect whether they have strong off-target effect. (boolean)
}
  \item{lower_bound_threshold}{
    lower bound on the average(or median) score of the seed to be considered effective (number)
}
  \item{higher_bound_threshold}{
    higher bound on the average(or median) score of the seed to be considered effective (number)
}
  \item{min_oligos_x_gene_threshold}{
    minimum number of siRNAs that each gene must have in order to be reported.
    Using this function we end up removing several seeds, and hence some genes remain with few oligos (1 or 2). You could think that this low number is not enough to be able to detect an effect on the gene, and hence you may like to remove these genes setting this variable appropriately.  (integer)
}
  \item{useMedian}{
    specify whether to use the mean or the median to compute the score of each seed (boolean)
}
  \item{removeGenes}{
	specify whether to remove the genes for which at least one siRNA targeting that gene has been found effective. This approach is suggested for the pooled libraries. (boolean)
}
  \item{include_current_gene}{
    lower bound on the average score of the seed to be considered effective (boolean)
}
  \item{progress_bar}{
    lower bound on the average score of the seed to be considered effective (boolean)
}


}
\value{
screen data frame where we removed the siRNAs that show a detectable seed effect. 
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 1000 rows
	uuk_screen_reduced = uuk_screen[1:1000,]

	screen_corrected = seed_removal(add_seed(uuk_screen_reduced))
}

