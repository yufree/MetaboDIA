\name{compare_ms2}
\alias{compare_ms2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compare two MS2 spectra
%%  ~~function to do ... ~~
}
\description{
Calculate a MS2 similarity score between two spectra.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compare_ms2(tmp, tmp.lib, ms2_tol, ms2_ppm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tmp}{
  First MS2 spectrum
%%     ~~Describe \code{tmp} here~~
}
  \item{tmp.lib}{
  Second MS2 spectrum
%%     ~~Describe \code{tmp.lib} here~~
}
  \item{ms2_tol}{
  MS2 mass tolerance in Dalton
%%     ~~Describe \code{ms2_ppm} here~~
}
  \item{ms2_ppm}{
  MS2 mass tolerance in ppm. ms2_tol will be ignored when this parameter is set by the user.
%%     ~~Describe \code{ms2_tol} here~~
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
\keyword{ MetaboDIA }
