\name{MetaProfiler-class}
\Rdversion{1.1}
\docType{class}
\alias{MetaProfiler-class}

\title{Class \code{"MetaProfiler"}}
\description{
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("MetaProfiler", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{design}:}{Object of class \code{"data.table"} ~~ }
    \item{\code{data}:}{Object of class \code{"data.table"} ~~ }
    \item{\code{incorporation_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{intensity_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{labeling_ratio_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{score_name}:}{Object of class \code{"character"} ~~ }
    \item{\code{peptide_column_no_PTMs}:}{Object of class \code{"character"} ~~ }
    \item{\code{peptide_column_PTMs}:}{Object of class \code{"character"} ~~ }
    \item{\code{protein_column}:}{Object of class \code{"character"} ~~ }
    \item{\code{peptide_centric}:}{Object of class \code{"logical"} ~~ }
    \item{\code{annotate_by_peptide}:}{Object of class \code{"character"} ~~ }
    \item{\code{master_tbl}:}{Object of class \code{"data.table"} ~~ }
    \item{\code{timepoints}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{time_zero}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{time_unit}:}{Object of class \code{"character"} ~~ }
    \item{\code{annotation_table}:}{Object of class \code{"data.table"} ~~ }
    \item{\code{taxon_columns}:}{Object of class \code{"character"} ~~ }
    \item{\code{function_columns}:}{Object of class \code{"character"} ~~ }
    \item{\code{deconvolve}:}{Object of class \code{"function"} ~~ }
    \item{\code{consts}:}{Object of class \code{"list"} ~~ }
    \item{\code{model}:}{Object of class \code{"function"} ~~ }
    \item{\code{pdf}:}{Object of class \code{"list"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "MetaProfiler" in the signature.
}
\references{
%%  ~~put references to the literature/web site here~~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%%  ~~objects to See Also as \code{\link{~~fun~~}}, ~~~
%%  ~~or \code{\linkS4class{CLASSNAME}} for links to other classes ~~~
}
\examples{
showClass("MetaProfiler")
}
\keyword{classes}
