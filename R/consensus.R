## This function builds a consensus library across the samples
##@param
#### ms2_spectrum: all the input MS2 spectra within CIU
#### ms2_tol: tolerance for ms2 fragment (in Dalton, it will be ignored if ms2_ppm is set)
#### ms2_ppm: tolerance for ms2 fragment (in ppm)
#### consus_filter: cutoff of the relative intensity of the consensus library
#### ms2_rep: cutoff of the reproducibility of the consensus peak across the samples in the library
#### ms2_IQR: inter-quartile-range of differences of the consensus peak and the corresponding peaks of each sample

consensus <- function(ms2_spectrum, ms2_tol, ms2_ppm, consus_filter, ms2_rep, ms2_IQR){
  n_spectrum <- length(ms2_spectrum)
  spectrum_all <- matrix(do.call(rbind,ms2_spectrum), ncol=2)
  spectrum_all <- matrix(spectrum_all[order(spectrum_all[ ,1]), ], ncol=2)

  ###Group ms2 spectra
  #difference between adjacent peaks
  if(nrow(spectrum_all)>1){
    mz_group <- list()
    mz_group[[1]] <- c(1)
    for(i in 2:nrow(spectrum_all)){
      mz_dif <- spectrum_all[i,1]-spectrum_all[(i-1),1]
      if(!is.null(ms2_ppm)) ms2_tol <- spectrum_all[(i-1),1] * ms2_ppm/ 1000000
      if(mz_dif < ms2_tol){
        mz_group[[length(mz_group)]] <- c(mz_group[[length(mz_group)]],i)
      } else {
        mz_group[[length(mz_group)+1]] <- i
      }
    }
    temp <- lapply(mz_group, function(x) {
      tmp <- matrix(spectrum_all[x,], ncol=2)
      if(length(x) > n_spectrum * ms2_rep){
        mz <- median(tmp[,1])
        mz_diff.i <- (tmp[,1] - mz) * 1000000 / mz
        if(IQR(mz_diff.i) <= ms2_IQR){
          inten <- mean(tmp[,2])
          list(lib.i=c(mz,inten), mz_diff.i=mz_diff.i)

        }
      }
    })
    id <- do.call("c", lapply(as.list(c(1:length(temp))), function(z) {
      temp2 <- temp[[z]]
      if(!is.null(temp2)){
        if(temp2$lib.i[2] >= consus_filter){
          return(z)
        }
      }
    }) )
    if(length(id) > 0){
      mz_var <- do.call(c,lapply(as.list(id), function(y) { temp[[y]]$mz_diff.i }))
      consensus_lib <- do.call(rbind,lapply(as.list(id), function(y){ temp[[y]]$lib.i }))
      consensus_lib[,2] <- consensus_lib[,2] / max(consensus_lib[,2])
      return(list(lib=matrix(consensus_lib, ncol=2), mz_diff=mz_var))
    } 
    else {
      return(list(lib=NULL, mz_diff=NULL))
    }

  } 
  else {
    return(list(lib=NULL, mz_diff=NULL))
  }

}
