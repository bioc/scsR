\name{transcribe_seqs}
\alias{transcribe_seqs}
\title{
transcribe_seqs
}
\description{
transcribe the sequences that are present in a specific column of the input data frame
}
\usage{
transcribe_seqs(df, seqColName="siRNA_seq", toDNA=FALSE, progress_bar=FALSE)
}
\arguments{
  \item{df}{
    input data frame.
}
  \item{seqColName}{
    name of the column that contains the sequences (character vector)
( each sequence must be in the format of a character vector, i.e. a simple string)
}
  \item{toDNA}{
    choose whether to transcribe to DNA (i.e. put T instead of U)  (boolean)
}
  \item{progress_bar}{
    choose whether disable printing warnings/messages. (boolean)
}
}
\value{
data frame given in input, but with the sequences transcribed.
}
\author{
Andrea Franceschini
}
\examples{
	data(uuk_screen)

	input_screen = head(uuk_screen, n=10)
	uuk_screen_transcribed = transcribe_seqs(input_screen)
}
