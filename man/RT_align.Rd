\name{RT_align}
\alias{RT_align}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Retention time alignment
}
\description{
Align CIUs across samples
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
RT_align(map_charge, RT_expand)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{map_charge}{
  CIU links across samples
%%     ~~Describe \code{map_charge} here~~
}
  \item{RT_expand}{
  Time window (in seconds) for aligning MS1 elution profiles across samples
%%     ~~Describe \code{RT_expand} here~~
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
