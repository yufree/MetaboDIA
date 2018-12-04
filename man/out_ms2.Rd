\name{out_ms2}
\alias{out_ms2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Output RT-aligned MS2 intensity table in the DIA-only framework
}
\description{
Output RT-aligned MS2 intensity table in the DIA-only framework
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
out_ms2(ms2_output, a, input, compound_list, f, ms2_tol, ms2_ppm, cv_n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms2_output}{
  Output file name
%%     ~~Describe \code{ms2_output} here~~
}
  \item{a}{
  Link between the formula and input data
%%     ~~Describe \code{a} here~~
}
  \item{input}{
  Input data information
%%     ~~Describe \code{input} here~~
}
  \item{compound_list}{
  Molecular formula list
%%     ~~Describe \code{compound_list} here~~
}
  \item{f}{
  Input file names
%%     ~~Describe \code{f} here~~
}
  \item{ms2_tol}{
  MS2 mass tolerance in Dalton
%%     ~~Describe \code{ms2_ppm} here~~
}
  \item{ms2_ppm}{
  MS2 mass tolerance in ppm. ms2_tol will be ignored when it is set
%%     ~~Describe \code{ms2_tol} here~~
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ MetaboDIA}
