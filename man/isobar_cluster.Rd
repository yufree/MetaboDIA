\name{isobar_cluster}
\alias{isobar_cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Cluster CIUs
}
\description{
Cluster the CIUs of identical formula at different RTs
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
isobar_cluster(file, n_pair, cut)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{file}{
  Output file name prior to clustering
%%     ~~Describe \code{file} here~~
}
  \item{n_pair}{
  Minimum number of pairs when calculating the correlation
%%     ~~Describe \code{n_pair} here~~
}
  \item{cut}{
  Average correlation cutoff between the clusters of CIUs
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
\keyword{ MetaboDIA }
