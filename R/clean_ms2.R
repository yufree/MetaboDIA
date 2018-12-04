## This function removes the possible noise in the MS2 spectrum
## @param
#### ms2: MS2 spectrum to be cleaned
#### min_basepeak: minimum intensity for basepeak, which is the peak with the highest intensity
#### min_int: minimum relative intensity scaled to the basepeak
#### min_noise: minimum intensity for noise signal

clean_ms2 <- function(ms2, min_basepeak, min_int, min_noise){
  if(nrow(ms2) > 1){
    noise <- min_noise
    basepeak <- max(ms2[,2])
    if(basepeak > min_basepeak) {
      noise <- max(noise, basepeak*min_int)
      ms2 <- matrix(ms2[which(ms2[,2]>noise), ], ncol=2)
      ms2[,2] <- ms2[,2] / basepeak
      return(list(basepeak=basepeak, ms2=ms2))
    } 
    else {
      return(list(basepeak=NULL, ms2=NULL))
    }
  } 
  else {
    return(list(basepeak=NULL, ms2=NULL))
  }
}
