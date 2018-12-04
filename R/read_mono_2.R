## This function reads the mono_peak_file, which is the output file from convert_cam_output function without
## considering of MS2 spectra (this is for MS1-only workflow)
##@param
#### mono_peak_file: file name
#### ms1_tol: tolerance for MS1 peak (in Dalton, it will be ignored if ms1_ppm is set)
#### ms1_ppm: tolerance for MS1 peak (in ppm)
#### intensity_type: peak apex height or area under the curve
read_mono_2 <- function(mono_peak_file, ms1_tol, ms1_ppm, intensity_type){
  ms1 <- read.csv(mono_peak_file)
  metab <- list()
  for(i in 1:nrow(ms1) ){
    mz_ms1 <- ms1$mz[i]

    RT <- ms1$rt[i]
    RT_range <- c(ms1$rtmin[i], ms1$rtmax[i])
    if(intensity_type=="area"){
      intensity_ms1 <- ms1$intb[i]
    } 
    else if(intensity_type=="height"){
      intensity_ms1 <- ms1$maxo[i]
    }
    charge <- NA
    ms2_info <- list()
    basepeak <- list()
    metab[[i]] <- list(mz_ms1=mz_ms1, RT=RT, RT_range=RT_range, intensity_ms1=intensity_ms1,
                       charge=charge, N_iso=1, DBcompound=NA, formula.all=NA,
                       ms2=ms2_info, basepeak=basepeak)
  }
  list(metab=metab)
}
