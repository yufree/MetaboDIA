%% Copyright (C) <2016>  Gengbo Chen <gengbo_chen@u.nus.edu> and Hyungwon Choi <hyung_won_choi@nuhs.edu.sg>, National University of Singapore.
%% This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or any later %%version.
%% This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%% See the GNU General Public %%License for more details.
%% You can obtain a copy of the GNU General Public License along with this program at <http://www.gnu.org/licenses/>
\name{MetaboDIA-package}
\alias{MetaboDIA-package}
\alias{MetaboDIA}
\docType{package}
\title{
A framework for spectral assay library construction and targeted quantitative data extraction in metabolomics with data independent acquisition
}
\description{
Data independent acquisition-mass spectrometry (DIA-MS) coupled with liquid chromatography is a promising approach for rapid, automatic sampling of MS/MS data in untargeted metabolomics. However, wide isolation windows in DIA-MS generate MS/MS spectra containing a mixed population of fragment ions together with their precursor ions. This precursor-fragment ion map in a comprehensive MS/MS spectral library is crucial for relative quantification of fragment ions uniquely representative of each precursor ion. However, existing reference libraries are not sufficient for this purpose since the fragmentation patterns of small molecules can vary in different instrument setups. Here we developed a bioinformatics workflow called MetaboDIA to build customized MS/MS spectral libraries using a user’s own data dependent acquisition (DDA) data and to perform MS/MS-based quantification with DIA data, thus complementing conventional MS1-based quantification. MetaboDIA also allows users to build a spectral library directly from DIA data in studies of a large sample size. 

}
\details{
\tabular{ll}{
Package: \tab MetaboDIA\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-06-23\cr
License: \tab Apache 2.0\cr
}
}
\author{
Gengbo Chen; Hyungwon Choi

Maintainer: Hyungwon Choi <hwchoi@nuhs.edu.sg>
}

\references{
}
\keyword{ package }
\keyword{MetaboDIA }
