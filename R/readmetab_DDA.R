## Unlike the function read_mono, this function reads all MS1 (including formula identification) 
## and MS2 information for those peaks with isotopes (one of the output files from convert_cam_output)
## in a sample and create a object for the sample
##@param
#### um_out: File name of the peaks with isotope peaks (output from convert_cam_output function)
#### db: The coresponding DB search result of the camera output.
#### mzxml_file: the mzXML file  **** change the variable name????
#### cutoff_file: FDR_cutoff file which was developped for metmatch (not used in this version)
#### ms1_tol: Mass tolerance for ms1 peak (in Dalton, it will be ignored if ms1_ppm is set)
#### ms1_ppm: Mass tolerance for ms1 peak (in ppm)
#### min_basepeak: Minimum intensity for basepeak, which is the highest intensity peak
#### min_int: Minimum relative intensity to the basepeak
#### min_noise: Minimum intensity for noise signal
#### min_spectrum: Minimum number of MS2 peaks in a spectrum (otherwise the MS2 epctrum will be discarded)
#### mode: Positive or negative mode
#### ms2_tol: Mass tolerance for MS2 fragment (in Dalton, it will be ignored if ms2_ppm is set)
#### ms2_ppm: Mass tolerance for MS2 fragment (in ppm)
#### intensity_type: Peak apex height or area under the curve
#### metmatch_FDR: Cutoff of metmatch FDR (developped for integration with metmatch, not used in this version)
#### best_hit: Use the best hit or all hit (developped for integration with metmatch, not used in this version)

readmetab_DDA <- function(um_out, db, mzxml_file, cutoff_file, ms1_tol, ms1_ppm, 
                        min_spectrum, min_basepeak, min_int, min_noise, mode, ms2_mz_tol,
                        metmatch_FDR, best_hit, intensity_type){

  PROTON_MASS = 1.007276
  ms2 <- read_mzxml(as.character(mzxml_file), min_spectrum)
  ms1_iden <- read.csv(db)
  ms1_iden <- ms1_iden[order(ms1_iden$index.Q),]
  ms1_iden <- ms1_iden[which(!is.na(ms1_iden$hit.DB)),]
  if(!is.na(metmatch_FDR)){
    tmp_cut <- read.table(cutoff_file,skip=1,header=F,sep=" ")
    co <- tmp_cut[which(as.numeric(gsub("alpha_","",tmp_cut[,1]))==metmatch_FDR),2]
    ms1_iden <- ms1_iden[which(ms1_iden$score.Q>=co),]
  }
  ms1 <- read.csv(um_out)
  ms1 <- ms1[order(ms1$id),]
  ms1id_iden <- NULL

  #ms2 matching
  #####
  ms2_match <- lapply(ms2, function(x){
    RT_ms2 <- x$RT
    ms1_mz <- x$ms1
    ms2_data <- x$ms2
    id <- which(RT_ms2<=ms1$rtmax & RT_ms2>=ms1$rtmin)
    if(length(id)>0){
      tmp <- ms1[id,]
      if(!is.null(ms1_ppm)) ms1_tol <- ms1_mz*ms1_ppm/1000000
      if(min(abs(tmp$mz-ms1_mz))<=ms1_tol){
        id <- which(abs(tmp$mz-ms1_mz)==min(abs(tmp$mz-ms1_mz)))
        id <- unique(tmp$id[id])
      } else {
        id <- NULL
      }
    }
    list(RT=RT_ms2, ms1=ms1_mz, id=id, ms2=ms2_data)
  })

  #rebuid  data structure

  metab <- list()
  for(i in 1:nrow(ms1_iden) ){
    tmpid <- which(ms1$id==ms1_iden$index.Q[i])
    mz_ms1 <- ms1$mz[tmpid[1]]

    RT <- ms1$rt[tmpid][1]
    RT_range <- c(ms1$rtmin[tmpid[1]],ms1$rtmax[tmpid[1]])
    if(intensity_type=="area"){
      intensity_ms1 <- ms1$intb[tmpid[1]]
    } else if(intensity_type=="height"){
      intensity_ms1 <- ms1$maxo[tmpid[1]]
    }
    charge <- ms1$charge[tmpid[1]]
    if(best_hit){
      DBcompound <- as.character(ms1_iden$best.hit[i])
    } else {
      DBcompound <- as.character(ms1_iden$hit.DB[i])
    }
    formula.all <- as.character(ms1_iden$hit.DB[i])

    n_isotope <- length(tmpid)
    ms2_info <- list()
    basepeak <- list()

    metab[[i]] <- list(mz_ms1=mz_ms1, RT=RT, RT_range=RT_range, intensity_ms1=intensity_ms1,
                     charge=charge, N_iso=n_isotope, DBcompound=DBcompound, formula.all=formula.all,
                     ms2=ms2_info, basepeak=basepeak)
  }

  ###linking MS2 and add MS2 to the structure record unmapped MS2

  with_ms2 <- NULL

  for(i in 1:length(ms2_match)) {
    if(length(ms2_match[[i]]$id)>0) {

      for(j in 1:length(ms2_match[[i]]$id)){
        ind <- which(ms1_iden$index.Q==ms2_match[[i]]$id[j])
        if(length(ind)>0){

          cl_m <- clean_ms2(ms2_match[[i]]$ms2,min_basepeak,min_int,min_noise)
          ms2_cleaned <- cl_m$ms2
          bp <- cl_m$basepeak

          if(length(ms2_cleaned)>0){
            metab[[ind]]$ms2[[length(metab[[ind]]$ms2)+1]] <- ms2_cleaned
            metab[[ind]]$basepeak[[length(metab[[ind]]$basepeak)+1]] <- bp
            with_ms2 <- c(with_ms2,ind)

          }
        }
      }
    }
  }
  ####output
  list(metab=metab, with_ms2=unique(with_ms2))
}
