\name{enrichment_heatmap}
\alias{enrichment_heatmap}
\title{
enrichment_heatmap
}
\description{
Produces an heatmap showing the enriched annotations that are found in the input vectors of gene identifiers.
}
\usage{
enrichment_heatmap(genesVectors, vectorsNames, output_file=NULL, title="", limit=400, species_ncbi_taxonomy_id=9606,
                              enrichmentType="Process", limitMultiPicture=NULL, fdr_threshold=0.05, pvalue_threshold=NULL, 
                              cexRow=NULL, cexCol=1, STRINGversion="9_05", selectTermsVector=NULL, iea = TRUE, sortingMethod="rowMeans")
}
\arguments{
  \item{genesVectors}{
    list containing several sorted vectors of genes (i.e. columns of the heatmap) for which to compute the enrichment in pathways (list)
}
  \item{vectorsNames}{
    names of the vectors (to be displayed as column labels on the heatmap) (vector of strings)
}
  \item{output_file}{
    path to an output file where to store the heatmap (this file should have the pdf extension) (character vector)
}
  \item{title}{
    title of the heatmap graph (character vector)
}
  \item{limit}{
    considers only the top genes in the vector (integer)
}
  \item{species_ncbi_taxonomy_id}{
    ncbi taxonomy id of the organism (e.g. 9606 for Human) (integer)
}
  \item{enrichmentType}{
    type of Enrichment of the heatmap (either Process or KEGG.  The first tests for enrichment in GO biological processes, while the second tests for the enrichment in KEGG pathways)  (character vector)
}
  \item{limitMultiPicture}{
    number of rows of the heatmap before to start a new page in the pdf (integer)
}
  \item{fdr_threshold}{
    considers only the rows with at least one element below this threshold (number)
}
  \item{pvalue_threshold}{
    considers only the rows with at least one element below this threshold (number)
}
  \item{cexRow}{
    size of the row labels (number)
}
  \item{cexCol}{
    size of the columns' labels (number)
}
 \item{STRINGversion}{
    specify the version of STRING to use for the enrichment annotations (by default 9_05) (character vector)
}
 \item{selectTermsVector}{
    specify the terms to select. Each term must fully contain at least one string of this vector.
	This parameter can be used when we want to limit the output of the method, for example to fit the output image in one page of an article (vector of strings).
}
 \item{iea}{
	specify whether to use Electronic Inferred Association annotations (to be used in case you are querying the GeneOntology). (boolean)
}
 \item{sortingMethod}{
 specify whether a sorting method should be applied. For the moment, the only available method is rowMeans.  (character vector) 
}
}
\value{
matrix that is used to generate the heat map 
}
\author{
Andrea Franceschini
}
\examples{
data(uuk_screen)
data(uuk_screen_dh)
 \dontrun{ 	
heatmapMatrix = enrichment_heatmap( list( uuk_screen$GeneID,
                          arrange(add_rank_col(uuk_screen), log_pval_rsa)$GeneID,
                          uuk_screen_dh$GeneID
),
                    list("Qiagen", "Qiagen (RSA)", "Dharmacon"),
                    limit=400,
                    enrichmentType = "Process",
                    output_file=NULL,
                    title="Uuk Cell Killers",
			selectTermsVector=c("cycle")
)

}
}