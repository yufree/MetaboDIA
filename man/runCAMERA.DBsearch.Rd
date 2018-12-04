\name{runCAMERA.DBsearch}
\alias{runCAMERA.DBsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run XCMS CAMERA and Database search (molecular formula assignment)
%%  ~~function to do ... ~~
}
\description{
Run XCMS CAMERA and Database search (molecular formula assignment)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runCAMERA.DBsearch(dir=getwd(), DB.file, adduct.file, mode, n_core=4, ppm=30, prefilter=c(2, 50), method="centWave", peakwidth=c(5, 60), cor_eic_th=0.6, scanrange=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dir}{
  Path to files
%%     ~~Describe \code{dir} here~~
}
  \item{DB.file}{
  Database file name
%%     ~~Describe \code{DB.file} here~~
}
  \item{mode}{
  Positive or negative ionization mode
%%     ~~Describe \code{mode} here~~
}
  \item{adduct.file}{
  Adduct ion information to be considered in molecular formula assignment
%%     ~~Describe \code{DDA_dir} here~~
}
  \item{n_core}{
  Number of cores for parallel computing
%%     ~~Describe \code{n_core} here~~
}
  \item{ppm}{
  ppm parameter used in xcmsSet
%%     ~~Describe \code{ppm} here~~
}
   \item{prefilter}{
  prefilter parameter used in xcmsSet
%%     ~~Describe \code{n.core} here~~
}
  \item{method}{
  peak detection method used in xcmsSet
%%     ~~Describe \code{n.core} here~~
}
  \item{peakwidth}{
  peakwidth parameter used in xcmsSet
%%     ~~Describe \code{n.core} here~~
}
  \item{cor_eic_th}{
  cor_eic_th parameter used in CAMERA
%%     ~~Describe \code{n.core} here~~
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
\keyword{MetaboDIA }
