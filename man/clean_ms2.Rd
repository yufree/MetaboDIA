\name{clean_ms2}
\alias{clean_ms2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Clean a MS2 spectrum
}
\description{
Convert the intensities to relative intensities.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clean_ms2(ms2,  min_basepeak, min_int, min_noise)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{ms2}{
  A table of m/z and intensity of MS2 peaks
%%     ~~Describe \code{ms2} here~~
}
  \item{ms1_tol}{
  MS1 mass tolerance in Dalton
%%     ~~Describe \code{ms1_tol} here~~
}
  \item{ms1_ppm}{
  MS1 mass tolerance in ppm. ms1_tol will be ignored when this parameter is set by the user.
%%     ~~Describe \code{ms1_ppm} here~~
}
  \item{min_basepeak}{
  Minimum intensity of basepeak
%%     ~~Describe \code{min_basepeak} here~~
}
  \item{min_int}{
  Minimum relative intensity of MS2 peak
%%     ~~Describe \code{min_int} here~~
}
  \item{min_noise}{
  Minimum noise level
%%     ~~Describe \code{min_noise} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{basepeak}{Intensity of basepeak}
  \item{ms2 }{m/z and relative intensity of MS2 data}
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
