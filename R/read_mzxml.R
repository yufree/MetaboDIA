## This function reads the MS2 information form mzXML file
##@param
####fname: file name
####min_spectrum: minimum number of MS2 peaks (otherwise the MS2 epctrum will be discarded)
read_mzxml <- function(fname, min_spectrum){
  library("mzR")
  OMSMSF <- openMSfile(fname)
  HOMSMSF <- header(OMSMSF)
  HOMSMSF1 <- HOMSMSF[HOMSMSF$msLevel==2,]
  aqn1 <- HOMSMSF1$acquisitionNum
  RT <- HOMSMSF1$retentionTime
  ms1 <- HOMSMSF1$precursorMZ
  ncharge <- HOMSMSF1$precursorCharge
  charge <- ncharge[ncharge==0]  <-  NA
  vals1 <- data.frame(aqn1=aqn1, RT=RT, ms1=ms1, charge=charge)
  MS2 <- lapply(1:length(aqn1), function(x){
    ele1 <- aqn1[x]
    PCL <- mzR::peaks(OMSMSF,ele1)
    if(nrow(PCL)>min_spectrum){
      list(RT = vals1$RT[x], ms1=vals1$ms1[x], charge=vals1$charge[x], ms2=PCL)
    } else {
      NULL
    }
  })
  ##filter the empty entries
  Filter(length, MS2)
}
