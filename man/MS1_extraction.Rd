\name{MS1_extraction}
\alias{MS1_extraction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extraction MS1 data with alignment
}
\description{
MS1 workflow: identify molecular formula to MS1 precursor ion peaks and align MS1 intensities across samples
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MS1_extraction(file_dir, n_core = 4, ms1_output = "MS1_extraction.txt", RT_expand = 6, cv_n=5, ms1_ppm=30,
 ms1_tol=0.001, metmatch_FDR=NA, best_hit=F, intensity_type="area")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file_dir}{
  Path to data
%%     ~~Describe \code{DDA_dir} here~~
}
  \item{n_core}{
  Number of cores used in parallel computing
%%     ~~Describe \code{n_core} here~~
}
  \item{ms1_output}{
  Output file name
%%     ~~Describe \code{ms1_output} here~~
}
  \item{RT_expand}{
  Time window (in seconds) for aligning MS1 elution profiles across samples
%%     ~~Describe \code{RT_expand} here~~
}
  \item{cv_n}{
  Data extraction: minimum number of non-missing extracted intensities
%%     ~~Describe \code{cut} here~~
}
  \item{ms1_tol}{
  MS1 mass tolerance in Dalton
%%     ~~Describe \code{ms1_ppm} here~~
}
  \item{ms1_ppm}{
  MS1 mass tolerance in ppm. ms1_tol will be ignored when it is set
%%     ~~Describe \code{ms1_tol} here~~
}
\item{metmatch_FDR}{
  The FDR cutoff if using MetMatch method. Set to "NA" if not using Metmatch method for compound identification.
%%     ~~Describe \code{DDA_dir} here~~
}
\item{best_hit}{
  Whether to use the best hit from MetMatch method. Default is "false"
  %%     ~~Describe \code{DDA_dir} here~~
}
\item{intensity_type}{
  Use "area" or "height" from XCMS output as the quantification measurements
  %%     ~~Describe \code{DDA_dir} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gengbo Chen; Hyungwon Choi
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  Please refer to the MetaboDIA user manual for detailed example

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MetaboDIA }
