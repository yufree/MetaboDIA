## This function outputs the MS1 intensities after building the consensus library
##@param
#### ms1_output: Output file name
#### lib: Consensus library(DDA or DIA based)
#### input: A data object containing all the processed information for samples
#### compound_list: Formula list of the consensus library
#### f: DIA file names
#### cv_n: The minimum number of samples to compute CV

out_ms1 <- function(ms1_output, lib, input, compound_list, f, cv_n){
  file <- ms1_output
  file1 <- sub(".txt", ".isotope.txt", file)
  fields <- c("formula", "charge", "RT","ref_mz", "CV", "missing", "mean_elution_time", as.vector(row.names(f)))
  write.table(t(fields),file, sep="\t", row.names=F, col.names=F, quote=F)
  write.table(t(fields),file1, sep="\t", row.names=F, col.names=F, quote=F)
  for(i in 1:length(lib)){
    tmp <- lib[[i]]
    for(j in 1:length(tmp)){
      tmp2 <- tmp[[j]]
      out <- c(tmp2$char,tmp2$rt)
      tmp_rts <- rep(NA, length(input))
      int <- rep(0, length(input))
      n_isoto <- rep(0, length(input))
      tmp_ref_mz<-rep(NA, length(input))
      for(k in 1:length(tmp2$alignment)){
        tmp3 <- tmp2$alignment[[k]]
        if(!is.null(tmp3)){
          isotope_counts <- sapply(tmp3[1,], function(x){
            input[[k]]$metab[[x]]$N_iso
          })
          elution_length <- tmp3[4,]-tmp3[3,]
          intensities <- sapply(tmp3[1,], function(x){
            input[[k]]$metab[[x]]$intensity_ms1
          })
          int[k] <- max(intensities)
          n_isoto[k] <- isotope_counts[which.max(intensities)]
          tmp_rts[k] <- elution_length[which.max(elution_length)]
          tmp_ref_mz[k]<-median(tmp3[6,])
        }
      }
      mean_elution <- mean(tmp_rts, na.rm=T)
      ref_mz<-median(tmp_ref_mz,na.rm = T)
      ####calculate CVs here
      tmp_cv <- log2(int[which(int>0)])
      if(length(tmp_cv)>=cv_n){
        cv_score <- sd(tmp_cv)/mean(tmp_cv)
      } else {
        cv_score <- NA
      }
      missing <- 1 - length(tmp_cv)/length(int)
      out <- c(out, ref_mz,cv_score, missing, mean_elution)
      write.table(t(c(out,int)), file, append=T, sep="\t", row.names=paste(compound_list[[i]],collapse = ";"), col.names=F, quote=F)
      write.table(t(c(out,n_isoto)), file1, append=T, sep="\t", row.names=paste(compound_list[[i]],collapse = ";"), col.names=F, quote=F)
    }
  }
}
