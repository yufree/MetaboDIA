## This function read the mono_peak_file, which is the output file from convert_cam_output function together
## with MS2 spectra, link the MS2 and MS1 information
##@param
#### mgf_file: MS2 file (mzXML file for DDA or mgf for DIA )
#### mono_peak_file: MS1 peaks without isotopes which one of the output files from convert_cam_output function
#### ms1_tol: Mass tolerance for ms1 peak (in Dalton, it will be ignored if ms1_ppm is set)
#### ms1_ppm: Mass tolerance for ms1 peak (in ppm)
#### intensity_type: Peak apex height or area under the curve
#### min_basepeak: Minimum intensity for basepeak, which is the highest intensity peak
#### min_int: Minimum relative intensity to the basepeak
#### min_noise: Minimum intensity for noise signal
#### min_spectrum: Minimum number of MS2 peaks (otherwise the MS2 epctrum will be discarded)
#### mode: Positive or negative ionization mode
#### ms2_tol: Mass tolerance for MS2 fragment (in Dalton, it will be ignored if ms2_ppm is set)
#### ms2_ppm: Mass tolerance for MS2 fragment (in ppm)
#### intensity_type: Peak apex height or area under the curve

read_mono <- function(mgf_file, mono_peak_file, ms1_tol, ms1_ppm, min_spectrum, min_basepeak,
                    min_int, min_noise, mode, ms2_mz_tol, intensity_type){
  PROTON_MASS = 1.007276

  ####if DDA file use readmzxml DIA usereadmgf
  if(length(mgf_file)==1){
    ms2 <- read_mzxml(mgf_file, min_spectrum)
  } else {
    ms2 <- c(read_mgf(mgf_file[1],min_spectrum), read_mgf(mgf_file[2],min_spectrum))
  }

  ms1 <- read.csv(mono_peak_file)

  #map with MS2 data
  ms2_match <- lapply(ms2,function(x){
    RT_ms2 <- x$RT
    ms1_mz <- x$ms1
    ms2_data <- x$ms2
    id <- which(RT_ms2<=ms1$rtmax & RT_ms2>=ms1$rtmin)
    if(length(id)>0){
      tmp <- ms1[id,]
      if(!is.null(ms1_ppm)) ms1_tol <- ms1_mz*ms1_ppm/1000000
      if(min(abs(tmp$mz-ms1_mz))<=ms1_tol){
        id2 <- which(abs(tmp$mz-ms1_mz)==min(abs(tmp$mz-ms1_mz)))
        id3 <- id[id2]
      } else {
        id3 <- NULL
      }
    } else {
      id3 <- NULL
    }
    list(RT=RT_ms2, ms1=ms1_mz, id=id3, ms2=ms2_data)
  })


  ####first read the MS1 information and creat a list for it
  metab <- list()
  for(i in 1:nrow(ms1) ){
    mz_ms1 <- ms1$mz[i]

    RT <- ms1$rt[i]
    RT_range <- c(ms1$rtmin[i],ms1$rtmax[i])
    if(intensity_type=="area"){
      intensity_ms1 <- ms1$intb[i]
    } else if(intensity_type=="height"){
      intensity_ms1 <- ms1$maxo[i]
    }
    charge <- NA
    ms2_info <- list()
    basepeak <- list()
    metab[[i]] <- list(mz_ms1=mz_ms1, RT=RT, RT_range=RT_range, intensity_ms1=intensity_ms1,
                       charge=charge, N_iso=1, DBcompound=NA, formula.all=NA, ms2=ms2_info, basepeak=basepeak)
  }

  ######identify the MS1 peaks with MS2 information, add MS2 information
  with_ms2 <- NULL
  for(i in 1:length(ms2_match)) {
    if(length(ms2_match[[i]]$id)>0) {
      for(j in 1:length(ms2_match[[i]]$id)){
        ind <- ms2_match[[i]]$id[j]
        cl_m <- clean_ms2(ms2_match[[i]]$ms2, min_basepeak, min_int, min_noise)
        ms2_cleaned <- cl_m$ms2
        bp <- cl_m$basepeak

        if(length(ms2_cleaned)>0){
          tmp_char <- ceiling((max(ms2_cleaned[,1]) - ms2_mz_tol)/metab[[ind]]$mz_ms1)
          tmp_char <- max(tmp_char,abs(metab[[ind]]$charge),na.rm = T)
          if(mode=="positive"){
            metab[[ind]]$charge <- tmp_char
          } else {
            metab[[ind]]$charge <-  -tmp_char
          }

          metab[[ind]]$ms2[[length(metab[[ind]]$ms2)+1]] <- ms2_cleaned
          metab[[ind]]$basepeak[[length(metab[[ind]]$basepeak)+1]] <- bp
          with_ms2 <- c(with_ms2,ind)

        }
      }
    }
  }

  #output
  list(metab=metab,with_ms2=unique(with_ms2))
}

