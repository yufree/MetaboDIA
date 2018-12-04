\name{consensus}
\alias{consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Wrapper function for building a consensus library.
}
\description{
Construct a consensus spectral library across multiple samples.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
consensus(ms2_spectrum, ms2_tol, ms2_ppm, consus_filter, ms2_rep, ms2_IQR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ms2_spectrum}{
  A list of MS2 spectra objects
%%     ~~Describe \code{ms2_spectrum} here~~
}
  \item{ms2_tol}{
  MS2 mass tolerance in Dalton
%%     ~~Describe \code{ms2_tol} here~~
}
  \item{ms2_ppm}{
  MS2 mass tolerance in ppm. ms2_tol will be ignored when this parameter is set by the user.
%%     ~~Describe \code{ms2_ppm} here~~
}
  \item{consus_filter}{
  Average relative intensity cutoff
%%     ~~Describe \code{consus_filter} here~~
}
  \item{ms2_rep}{
  Cross-sample reproducibility threshold for a specific MS2 peak in the consensus library
%%     ~~Describe \code{ms2_rep} here~~
}
  \item{ms2_IQR}{
  Inter-quartile range of the peaks aligned to the consensus peak (in ppm)

%%     ~~Describe \code{ms2_IQR} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{lib }{Consensus library object}
  \item{mz_diff }{A vector of m/z differences between and MS2 peaks and the consensus peaks which they are aligned to.}
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
\keyword{ MetaboDIA }
