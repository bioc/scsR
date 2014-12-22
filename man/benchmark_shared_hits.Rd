\name{benchmark_shared_hits}
\alias{benchmark_shared_hits}
\title{
benchmark_shared_hits
}
\description{
This method can be used to benchmark sorted gene vectors (A) that comes out from a siRNA screen.
The benchmark is done against other sorted gene vectors (B) that we know to contain high density of real hits (e.g. the results of a second siRNA screen performed with a different library).
The benchmark is performed simply comparing the top n hits of the two lists. If the two lists contain many shared best hits than we have a strong statistical signal.
Then we display the number of shared best hits for different n, in a graph (if visualize_pval variable is set to true the pvalue of the t-test is plotted instead of the number of shared hits).
}
\usage{
benchmark_shared_hits(glA, glB, col, avoidIntersectL=FALSE, 
                                     output_file=NULL, npoints=400, title="", scaleAXPoint = 1, 
                                     scaleBXPoint = NULL, fixedBXPoint=400, displayRandomMultipleLines=TRUE, 
                                     nrandom=20, intersectGenes=TRUE, visualize_pval=FALSE, max_ylim=NULL, xlab=NULL, ylab="shared hits")
}
\arguments{
  \item{glA}{
    sorted list containing one or more sorted vectors of genes (i.e. hits of a genome wide screen sorted by significance).
    Each element i of glA will be benchmarked against element i of glB.  In case glB contains only one element, each glA vector will be benchmarked against glB[1].
}
\item{glB}{
    sorted list containing one or more sorted vectors of genes (i.e. hits of a genome wide screen sorted by significance).
}
  \item{col}{
    sorted vector of booleans (a boolean i in the vector corresponds to the shared hits of glA[i] with glB[i] )
    
}
  \item{avoidIntersectL}{
    sorted vector of colors (a color i in the vector corresponds to the shared hits line obtain intersecting glA[i] with glB[i] )
    To perform the benchmark we construct a background to be used (this background is given by the intersection of all the glA and glB vectors)
    When an element i of the vector is set to TRUE, we don't use the elements of glA[i] to compute the vector.
    This allows to benchmark also methods that do output only few putative good genes (instead of a sorted list of all the genes tested).
}
  \item{npoints}{
    number of points on the x-axis of the graph (integer)
}
  \item{nrandom}{
    number of random lines to compute (in order to infer the variation of the noise) (integer)
}
  \item{output_file}{
    path to the output file where to store the graph (character vector)
}
\item{title}{
    title of the graph (character vector)
}
\item{scaleAXPoint}{
   for position x in the graph we compare the best x * scaleAXPoint best hits of the genesA vector (integer)
}
\item{scaleBXPoint}{
    for position x in the graph we compare the best x * scaleBXPoint best hits of the genesB vector (integer)
}
\item{fixedBXPoint}{
    for position x in the graph we compare the best fixedBXPoint best hits of the genesB vector (integer)
}
\item{intersectGenes}{
    specify whether to intersect the genes from the various input vectors to form a suitable background to be used for the benchmark. (boolean)
}
\item{visualize_pval}{
    specify whether a p-value (derived by an hypergeometric test) should be visualized instead of the number of shared hits. (boolean)
}
\item{displayRandomMultipleLines}{
    specify whether to display several random lines in the graph (instead of only one line that is the average of all the random lines) (boolean)
}
\item{max_ylim}{
    y upper limit (integer)
}
\item{xlab}{
    xlab (character vector)
}
\item{ylab}{
    ylab (character vector)
}

}
\author{
Andrea Franceschini
}
\examples{
data(uuk_screen)
data(uuk_screen_dh)

benchmark_shared_hits(
  glA=list(
    uuk_screen[1:1000,]$GeneID, 
    arrange(add_rank_col(uuk_screen[1:1000,]), log_pval_rsa)$GeneID
  ),
  glB=list(uuk_screen_dh$GeneID),
  col=c("black", "blue"),
  title="UUKUNIEMI Hela Cell Killers"
)
}
