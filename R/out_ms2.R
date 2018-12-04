## This function outputs the MS2 intensity data for DIA_DIA_workflow
##@param
#### ms2_output: Output file name
#### lib: DIA based consensus ilbrary
#### input: A data object that contain all the processed information for DIA samples
#### f: DIA file names
#### ms2_tol: Mass tolerance for MS2 fragment (in Dalton, it will be ignored if ms2_ppm is set)
#### ms2_ppm: Mass tolerance for MS2 fragment (in ppm)
#### cv_n: The minimum number of samples to compute CV

out_ms2 <- function(ms2_output, lib, input, compound_list, f, ms2_tol, ms2_ppm, cv_n){
  file <- ms2_output
  ####output header
  fields <- c("formula", "charge", "RT", "ms2_lib_mz", "ms2_lib_intensity", "CV", "missing", as.vector(row.names(f)))
  write.table(t(fields),file, sep="\t", row.names=F, col.names=F, quote=F)
  #########output each alignment
  for(i in 1:length(lib)){
    ###each compound
    tmp <- lib[[i]]
    for(j in 1:length(tmp)){
      #########each ms1 alignment
      tmp2 <- tmp[[j]]
      if(!is.null(tmp2[[2]])) {
        con_lib <- tmp2$consensus  #####if there is a consensus library
        ms2_align <- matrix(0,nrow=nrow(con_lib), ncol=length(input))
        charge <- tmp2$char
        rt <- tmp2$rt

        for(k in 1:length(tmp2$alignment)){
          tmp3 <- tmp2$alignment[[k]]
          ######there is ms2 in sample k
          if(!is.null(tmp3)){
            ms2 <- NULL
            ms2_id <- intersect(tmp3[1,], input[[k]]$with_ms2)
            if(length(ms2_id)>1){
              ms2 <- do.call("rbind",lapply(as.list(ms2_id), function(x){
                ms2 <- input[[k]]$metab[[x]]$ms2
                if(length(ms2)>1) {
                  ms2_combine <- NULL
                  for(w in 1:length(ms2)){
                    ms2_tmp=ms2[[w]]
                    ms2_bp <- input[[k]]$metab[[x]]$basepeak[[w]]
                    ms2_tmp[,2] <- ms2_tmp[,2]*ms2_bp
                    ms2_combine <- rbind(ms2_combine,ms2_tmp)
                  }
                  return(ms2_combine)
                } else {
                  ms2=ms2[[1]]
                  ms2_bp <- input[[k]]$metab[[x]]$basepeak[[1]]
                  ms2[,2] <- ms2[,2]*ms2_bp
                  return(ms2)
                }
              }))
            } else if (length(ms2_id)==1){
              ms2 <- input[[k]]$metab[[ms2_id]]$ms2
              if(length(ms2)>1) {
                ms2_combine <- NULL
                for(w in 1:length(ms2)){
                  ms2_tmp=ms2[[w]]
                  ms2_bp <- input[[k]]$metab[[ms2_id]]$basepeak[[w]]
                  ms2_tmp[,2] <- ms2_tmp[,2]*ms2_bp
                  ms2_combine <- rbind(ms2_combine,ms2_tmp)
                }
                ms2 <- ms2_combine
              } else {
                ms2=ms2[[1]]
                ms2_bp <- input[[k]]$metab[[ms2_id]]$basepeak[[1]]
                ms2[,2] <- ms2[,2]*ms2_bp
              }
            }
            ######align ms2 with con_lib if there is ms2
            if (length(ms2_id)>0){
              ms2_align[,k] <- sapply(con_lib[,1], function(x){
                if(!is.null(ms2_ppm)) ms2_tol=x*ms2_ppm/1000000
                id <- which(abs(ms2[,1]-x)<=ms2_tol)
                if(length(id)==0){
                  return (0)
                } else if (length(id)==1){
                  return (ms2[id,2])
                } else {
                  return (max(ms2[id,2]))
                }
              })
            }

          }
        }
        #############ms2 extraction done and start to check cv and output
        cv_score <- apply(ms2_align,1,function(xx){
          tmp_cv <- log2(xx[which(xx>0)])
          if(length(tmp_cv) >= cv_n){
            sd(tmp_cv)/mean(tmp_cv)
          } else {
            return (NA)
          }
        })
        missing <- apply(ms2_align,1,function(xx){
          length(which(xx>0))/length(xx)
        })
        out <- cbind(charge, rt, matrix(con_lib,nrow = length(cv_score)), 
                     cv_score, missing, matrix(ms2_align,nrow = length(cv_score)))
        write.table(out, file, append = T, sep="\t", 
                    row.names=rep(paste(tmp2$formula[,1],collapse = ";"), nrow(out)), col.names=F, quote=F)
      }
    }
  }
}
