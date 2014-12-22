\name{plot_seeds_methods}

\alias{plot_effective_seeds_head}
\alias{plot_seeds_oligo_count}
\alias{plot_seed_score_sd}
\alias{plot_screen_seeds_count}

\title{
  plot seeds utility methods
}

\description{
  Plots informations about the effect of the seed on the screen
}

\details{
  The methods currently depend on the type of \code{to}:
  \describe{
    \item{plot_effective_seeds_head}{
	barplot that represents the most effective seeds as bar (the length of the bars corresponds to the average score of the oligos that contain that seed). A background bar is shown under every seed. We obtain these bar simply randomizing the score column of the screen (and they well represent the noise level).
	}
    \item{plot_seeds_oligo_count}{
For each seed that is found in the siRNA screen, plots the number of oligos that contain that seed.
    }
\item{plot_seed_score_sd}{
For each seed plot its average score and its standard deviation.
}
\item{plot_screen_seeds_count}{
For each siRNA oligo, plot the number of the other oligos in the screen that share the same seed.
}
  }
}


\usage{
plot_effective_seeds_head(screen, seedColName="seed7",  scoreColName="score", enhancer_analysis=FALSE, min_oligos_x_seed=10, number_of_seeds=20, output_file=NULL, color="#CCCCCC33", colorBG="#0000CC11", xlim=c(-4,4), title="")

plot_seeds_oligo_count(screen, seedColName="seed7", scoreColName="score", output_file=NULL)
}

\arguments{
  \item{screen}{
    data frame containing the results of the siRNA experiment.
}
  \item{seedColName}{
    name of the column that contains the seed of the siRNA oligo sequences of the screen (character vector)
(the sequences have to be provided in the guide/antisense orientation and each sequence must be in the format of a character vector, i.e. a simple string)
}
  \item{scoreColName}{
    name of the column that contains the score of the screen (character vector)
}
  \item{enhancer_analysis}{
    if set to true plot the seeds that cause the oligos to have an higher score instead of a lower score. (boolean)
}
\item{min_oligos_x_seed}{
    minimum number of oligos that seed must have in order to be considered (integer)
}
\item{number_of_seeds}{
    maximum number of seeds to represent in the graph (by default the top 20 seeds are shown) (integer)
}
\item{output_file}{
    name of the pdf file where to store the graph (character vector)
}
\item{color}{
    color of the bars that represent the seeds (character vector)
}
\item{colorBG}{
    color of the bars that represent the noise (i.e. analysis executed on randomized data) (character vector)
}
\item{xlim}{
    xlim of the graph (number)
}
\item{title}{
    title of the graph (number)
}
}

\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	# to speed up the example we use only the first 2500 rows
	uuk_screen_reduced = uuk_screen[1:5000,]

	plot_effective_seeds_head(add_seed(uuk_screen_reduced))
}
