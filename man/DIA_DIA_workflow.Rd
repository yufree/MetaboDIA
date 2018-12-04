\name{DIA_DIA_workflow}
\alias{DIA_DIA_workflow}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Alternative workflow in MetaboDIA
}
\description{
DIA-DIA framework:align the MS intensities of DIA data. generate consensus library and extract the MS2 intensities from DIA data

%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DIA_DIA_workflow(DIA_dir, n_core=4, mode, DB.file, adduct.file, min_sample=0.2, consensus_filter=0.05, ms2_rep=0.5, ms2_IQR=80, RT_expand=6, ms1_ppm=30, ms1_tol=0.001, ms2_ppm=40, ms2_tol=0.01, min_noise=20, min_int=0.01, min_basepeak=100, min_spectrum=1, ms2_output="DIA_ms2.txt", r_score=0.5, cv_n=5, de_ppm=40, ms2_mz_tol, metmatch_FDR=NA, best_hit=F, intensity_type="area")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{DIA_dir}{
  Path to DIA data
%%     ~~Describe \code{DIA_dir} here~~
}
\item{n_core}{
  Number of cores used in parallel computing
%%     ~~Describe \code{n_core} here~~
}
\item{mode}{
  "positive" or "negative"
%%     ~~Describe \code{DDA_dir} here~~
}
  \item{DB.file}{
  Compound identification database
%%     ~~Describe \code{DDA_dir} here~~
}
  \item{adduct.file}{
  Adduct ion information to be considered in molecular formula assignment
%%     ~~Describe \code{DDA_dir} here~~
}

  \item{identification.ppm}{
  The mass tolerance in ppm when doing the additional formula identification
%%     ~~Describe \code{DDA_dir} here~~
}
  \item{min_sample}{
  Library construction: minimun fraction of samples that have MS2 data in order to build a MS2 consensus library
%%     ~~Describe \code{min_sample} here~~
}
  \item{consensus_filter}{
  Library construction: the average relative intensity cutoff
%%     ~~Describe \code{consensus_filter} here~~
}
  \item{ms2_rep}{
  Library construction: the cutoff for reproducibility across the MS2 spectra for a specific ms2 peak in the consensus library
%%     ~~Describe \code{ms2_rep} here~~
}
  \item{ms2_IQR}{
  Library construction: cutoff for the inter-quartile range of the peaks aligned to the consensus peak in ppm

%%     ~~Describe \code{ms2_IQR} here~~
}
  \item{RT_expand}{
  RT alignment: time window (in seconds) for aligning MS1 elution profiles across samples
%%     ~~Describe \code{RT_expand} here~~
}
  \item{ms1_tol}{
  MS1 mass tolerance in Dalton
%%     ~~Describe \code{ms1_ppm} here~~
}
  \item{ms1_ppm}{
  MS1 mass tolerance in ppm. ms1_tol will be ignored when this parameter is set by the user.
%%     ~~Describe \code{ms1_tol} here~~
}
  \item{ms2_tol}{
  MS2 mass tolerance in Dalton
%%     ~~Describe \code{ms2_ppm} here~~
}
  \item{ms2_ppm}{
  MS2 mass tolerance in ppm. ms2_tol will be ignored when this parameter is set by the user.
%%     ~~Describe \code{ms2_tol} here~~
}
  \item{min_noise}{
  Minimum noise level. Intensities below this value are ignored
%%     ~~Describe \code{min_noise} here~~
}
  \item{min_int}{
  Minimum relative intensity of MS2 peak
%%     ~~Describe \code{min_int} here~~
}
  \item{min_basepeak}{
  Minimum intensity of basepeak
%%     ~~Describe \code{min_basepeak} here~~
}
  \item{min_spectrum}{
  Minimum number of peaks in MS2 spectum
%%     ~~Describe \code{min_spectrum} here~~
}
  \item{RT_tol}{
  DIA data Extraction: Retention time tolerance (in seconds) when extracting DIA MS2 with a DDA library
%%     ~~Describe \code{RT_tol} here~~
}
  \item{ms2_output}{
  Output file name for MS2 aligned intensities
%%     ~~Describe \code{ms2_output} here~~
}
  \item{r_score}{
  Consensus library building: minimum similarity score between the consensus library and each spectrum
%%     ~~Describe \code{cut} here~~
}
  \item{cv_n}{
  Data extraction: minimum number of non-missing extracted intensities
%%     ~~Describe \code{cut} here~~
}
  \item{de_ppm}{
  Library construction: m/z tolerance in ppm when de-isotoping
%%     ~~Describe \code{cut} here~~
}
\item{ms2_mz_tol}{
  m/z tolerance when cleaning a ms2 spectrum
%%     ~~Describe \code{DDA_dir} here~~
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

}
\keyword{ MetaboDIA}
