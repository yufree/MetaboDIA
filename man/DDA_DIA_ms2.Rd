\name{DDA_DIA_ms2}
\alias{DDA_DIA_ms2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Reporting MS2 quantification results from the DDA-DIA pipeline
}
\description{
Output extracted MS2 quantification results
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DDA_DIA_ms2(file, DDA_lib, input, DDA_DIA_compound_map, com_map, compound_list_DDA, f, RT_expand, ms2_tol, ms2_ppm,cv_n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
  Output file name
%%     ~~Describe \code{file} here~~
}
  \item{DDA_lib}{
  DDA consensus library
%%     ~~Describe \code{DDA_lib} here~~
}
  \item{input}{
  DIA data
%%     ~~Describe \code{input} here~~
}
  \item{DDA_DIA_compound_map}{
  Link between the DDA formulas to DIA data
%%     ~~Describe \code{DDA_DIA_compound_map} here~~
}
  \item{com_map}{
  link between DIA formulas to DIA data
%%     ~~Describe \code{com_map} here~~
}
  \item{compound_list_DDA}{
  Formula list in DDA library
%%     ~~Describe \code{compound_list_DDA} here~~
}
  \item{f}{
  DIA input files
%%     ~~Describe \code{f} here~~
}
  \item{RT_expand}{
  Time window (in seconds) for aligning MS1 elution profiles across samples
%%     ~~Describe \code{RT_expand} here~~
}
 \item{ms2_ppm}{
  MS2 mass tolerance in Dalton
%%     ~~Describe \code{ms2_ppm} here~~
}
  \item{ms2_tol}{
  MS2 mass tolerance in ppm. ms2_tol will be ignored when this parameter is set by the user.
%%     ~~Describe \code{ms2_tol} here~~
}
  \item{RT_tol}{
  Retention Time tollance (seconds) when extracting DIA MS2 with DDA library
%%     ~~Describe \code{RT_tol} here~~
}
  \item{cv_n}{
  Data extraction: minimum number of non-missing extracted intensities
%%     ~~Describe \code{cut} here~~
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
}
\keyword{ MetaboDIA }
