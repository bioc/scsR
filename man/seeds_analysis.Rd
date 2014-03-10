\name{seeds_analysis}
\alias{seeds_analysis}
\title{
seeds_analysis
}
\description{
Create a data frame with several information on the effect of each seed in the genome-wide siRNA screen.
The average score of that seed is reported, together with the number of oligos that contain that seed.
Besides, suitable statistics are performed in order to estimate the p-value that the seed has an effect on the phenotype:
- Hypergeometric test (i.e. the probability that the seeds has more hits than expected by chance)
- Kolmogorov Smirnov test (i.e. the probability that you can obtain such high oligo scores by chance sampling from the entire score vector in the screen).
In addition we also report the human miRNAs that have the same seeds as the oligos (given that you could be interested to test them in the lab).
}
\usage{
seeds_analysis(screen, seedColName="seed7",  scoreColName="score", hit_th_val=NULL, 
                               enhancer_analysis=FALSE, spAvgColName=NULL, 
                               minCount=NULL, ks_enabled=FALSE, miRBase=NULL)
}
\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{seedColName}{
    name of the column that contains the seeds of the siRNA oligo sequences. (character vector)
(the sequences have to be provided in the guide/antisense orientation and each sequence must be in the format of a character vector, i.e. a simple string)
}
  \item{scoreColName}{
    name of the column that contains the score of the screen (character vector)
}
  \item{hit_th_val}{
    if the score of an oligo is below this threshold we define it as an hit.  
    This is then used to compute a p-value with an hypergeometric test for the seed.
    If this value is left to NULL, the best 10 percent of the oligos are considered as hits.
	(number)
}
  \item{enhancer_analysis}{
    specify the direction of the analysis.
    When this variable is set to FALSE it means that we are looking at the seeds that decrease the score of the oligos (when it is set to TRUE, it means we are looking at the seeds that increase the score of the oligos). 
	(booleam)
}
  \item{spAvgColName}{
    it is possible to specify the name of one column on which we want to perform the average, when we group for the seeds (for example, other than looking at the phenotype we may want to know the effect on the seed also on the cell number and display it on the same table). (vector of strings)
}
  \item{minCount}{
    minimum number of oligos in which a seed must be present in order to be reported in the output table (integer)
}
  \item{ks_enabled}{
    specify whether you want to compute also a Kolmogorov Smirnov test on the score of the seed. (boolean)
}
  \item{miRBase}{
    data frame object containing the human miRNAs and their sequences. The names of the columns must be "miRNA" and "seq" (data frame)
}
}
\value{
return a data frame with several information on the effect of each seed in the genome-wide siRNA screen.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)
	data(miRBase_20)

	# to speed up the example we use only the first 1000 rows
	uuk_screen_reduced = uuk_screen[1:1000,]

	seeds = seeds_analysis( add_seed(uuk_screen_reduced), miRBase = miRBase_20)
}