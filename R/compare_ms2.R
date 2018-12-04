## This function calculates a similarity score between two MS2 spectra 
##@params
#### tmp: MS2 spectrum 1
#### tmp.lib: MS2 spectrum 2
#### ms2_tol: tolerance for ms2 fragment (in Dalton, it will be ignored if ms2_ppm is set)
#### ms2_ppm: tolerance for ms2 fragment (in ppm)

compare_ms2<-function(tmp, tmp.lib, ms2_tol, ms2_ppm){
  ms2 <- tmp
  ms2[,2] <- ms2[,2] / max(ms2[,2])

  lib.ms2 <- tmp.lib
  lib.ms2[,2] <- lib.ms2[,2] / max(lib.ms2[,2])
  de <- sqrt(sqrt(sum(ms2[,2]^2)) * sqrt(sum(lib.ms2[,2]^2))) * max(sum(ms2[,2]), sum(lib.ms2[,2]))

  temp <- lapply(as.list(1:nrow(ms2)), function(y) {
    mz.dif <- abs(ms2[y,1] - lib.ms2[,1])
    if(!is.null(ms2_ppm)) ms2_tol = ms2[y,1] * ms2_ppm / 1000000
    if(min(mz.dif) > ms2_tol) {
      return (c(NA,NA,NA))
    } else {
      id <- which.min(mz.dif)
      delt.mass <- (ms2[y,1] - lib.ms2[id,1]) / ms2[y,1] * 1000000
      return(c(delt.mass, ms2[y,2], lib.ms2[id,2]))
    }
  })
  tmp2 <- do.call(rbind,temp)
  nu <- sqrt(sum(tmp2[,2] * tmp2[,3], na.rm = T)) * (sum(tmp2[,2], na.rm=T) + sum(tmp2[,3], na.rm=T)) / 2
  nu / de
}

