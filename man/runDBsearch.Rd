\name{runDBsearch}
\alias{runDBsearch}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run Database search (molecular formula assignment)
%%  ~~function to do ... ~~
}
\description{
Run Database search (molecular formula assignment)
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runDBsearch(dir=getwd(), DB.file, adduct.file, mode, n_core=4, ppm=30)
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
  \item{adduct.file}{
  Adduct ion information to be considered in molecular formula assignment
%%     ~~Describe \code{DDA_dir} here~~
}
  \item{mode}{
  positive or negative ioinization mode
%%     ~~Describe \code{mode} here~~
}
  \item{decoy}{
  Whether to use the decoy database in MetMatch
%%     ~~Describe \code{decoy} here~~
}
  \item{n_core}{
  Number of cores for parallel computing
%%     ~~Describe \code{n_core} here~~
}
  \item{ppm}{
  ppm parameter used in xcmsSet
%%     ~~Describe \code{ppm} here~~
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
