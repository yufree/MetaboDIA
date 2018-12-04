## Similar to readmetab_DDA function, this function does not match 
## MS2 information. It is developped for MS1_workflow only.
##@param
#### um_out: File name of the peaks with isotope peaks (output from convert_cam_output function)
#### db: The coresponding DB search result of the camera output.
#### cutoff_file: FDR_cutoff file which was developped for metmatch (not used in this version)
#### intensity_type: Peak apex height or area under the curve
#### metmatch_FDR: Cutoff of metmatch FDR (developped for integration with metmatch, not used in this version)
#### best_hit: Use the best hit or all hit (developped for integration with metmatch, not used in this version)

readmetab_DDA_2 <- function(um_out, db, cutoff_file, metmatch_FDR, best_hit, intensity_type){
  ms1_iden <- read.csv(db)
  ms1_iden <- ms1_iden[order(ms1_iden$index.Q), ]
  ms1_iden <- ms1_iden[which(!is.na(ms1_iden$hit.DB)), ]
  if(!is.na(metmatch_FDR)){
    tmp_cut <- read.table(cutoff_file, skip=1, header=F, sep=" ")
    co <- tmp_cut[which(as.numeric(gsub("alpha_","",tmp_cut[,1])) == metmatch_FDR), 2]
    ms1_iden <- ms1_iden[which(ms1_iden$score.Q >= co), ]
  }
  ms1 <- read.csv(um_out)
  ms1 <- ms1[order(ms1$id),]
  ms1id_iden <- NULL

  #rebuid  data structure
  metab <- list()
  for(i in 1:nrow(ms1_iden) ){
    tmpid <- which(ms1$id == ms1_iden$index.Q[i])
    mz_ms1 <- ms1$mz[tmpid[1]]

    RT <- ms1$rt[tmpid][1]
    RT_range <- c(ms1$rtmin[tmpid[1]], ms1$rtmax[tmpid[1]])
    if(intensity_type == "area"){
      intensity_ms1 <- ms1$intb[tmpid[1]]
    } 
    else if(intensity_type == "height"){
      intensity_ms1 <- ms1$maxo[tmpid[1]]
    }
    
    charge <- ms1$charge[tmpid[1]]
    if(best_hit){
      DBcompound <- as.character(ms1_iden$best.hit[i])
    } 
    else {
      DBcompound <- as.character(ms1_iden$hit.DB[i])
    }
    formula.all <- as.character(ms1_iden$hit.DB[i])

    n_isotope <- length(tmpid)
    ms2_info <- list()
    basepeak <- list()
    metab[[i]] <- list(mz_ms1=mz_ms1, RT=RT, RT_range=RT_range, 
                      intensity_ms1=intensity_ms1, charge=charge,
                      N_iso=n_isotope, DBcompound=DBcompound, 
                      formula.all=formula.all, ms2=ms2_info, basepeak=basepeak)
  }

  ####output
  list(metab=metab)
}
