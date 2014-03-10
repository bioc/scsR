\name{seed_correction}
\alias{seed_correction}
\title{
seed_correction
}
\description{
This method assumes that the seed effect acts in an "additive" way to the on-target signal. For example if we have an oligo X that has score -2, the method computes the average score of the other oligos in the libraries that contain the same seed as the oligo X. If for example this average score turns out to be -1.5, we can just subtract this score to the original oligo score to obtain the new "corrected" score ( (-2) - (-1.5) = -0.5 ).
However, the method assumes also that the correction factor ( -1.5 in the previous example ) should be multiplied by a coefficient c that reflects "how much" we suppose the effect is really additive  ( (-2) - c*(-1.5) ).  
The coefficient c can be a constant (e.g. 0.5 ) or it can vary depending on the behavior of the oligos that share the seed. In particular, we observed the the last approach to be the most successful.  Hence for our algorithm we set c = 0.4 + s.
s is a factor proportional to the distance of the standard deviation of the oligos that share the same seed with respect to the standard deviation that is expected, given their average score. This is because we observed that the expected standard deviation of the oligos that share a seed strictly depends on the average score as it can be seen using our plot-seed-score-sd function.
In particular s = sd_correction_coeff * quantile_std (sd_correction_coeff is a constant, by default set to 0.6, and quantile_std is the quantile of the standard deviation of the seeds that have an average score in the same interval as that of the oligos having the seed of the oligo X ).


}
\usage{
seed_correction(screen, seedColName="seed7", scoreColName="score", 
                                geneColName="GeneID", fixed_correction_coeff=0.4, 
                               sd_correction_coeff=0.6, min_siRNAs_x_seed=3, progress_bar=FALSE)
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
    name of the column that contains the scores of the screen. (character vector)
}
  \item{geneColName}{
    name of the column that contains the identifiers of the genes in the screen (character vector)
}
  \item{fixed_correction_coeff}{
      This coefficient is summed to the sd_correction_coefficient to obtain the final coefficient that is multiplyed to the correction factor (this final number must always be between 0 and 1 ).  (character vector)
}
  \item{sd_correction_coeff}{
    correction coefficient that is calibrated using the standard deviation of the oligos (i.e. oligos having a much lower standard deviation than expected are corrected multiplying the correction factor exactly by this coefficient... while oligos that have a normal or high standard deviation are corrected multiplying the correction factor by a number that is always lower than this variable, and that depends to the quantile of the standard deviation).
    This coefficient (after calibration) is summed to the fixed coefficient to obtain the final coefficient that is multiplyed to the correction factor (this final number must always be between 0 and 1 ).  (number)
}
  \item{min_siRNAs_x_seed}{
    This variable specify the minimum number of oligos that need to be present in the screen with the same seed as the oligo that we are correcting (the oligo that we are correcting is not included in the count).  (number)
}
  \item{progress_bar}{
    set this parameter to TRUE to show a progress bar. (boolean)
}
}
\value{
screen data frame with the score of the oligos corrected for the seed effect.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# To reduce the execution time in the example we trim the real dataset to contain only the first 500 rows
	# However in any real case the entire content of a genome-wide screen should be provided in as input.
	screen=uuk_screen[1:500,]

	screen_corrected = seed_correction(add_seed(screen))
}
