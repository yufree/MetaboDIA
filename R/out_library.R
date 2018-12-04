## This function outputs the consensus library to a file
##@param
#### lib: consensus library
#### file: file name

out_library <- function(lib, file, DB){
  fname <- file
  qc_score_out <- NULL
  write.table(NULL, fname, row.names=F, col.names=F, quote=F)
  for (i in 1:length(lib)){

    tmp <- lib[[i]]

    for(j in 1:length(tmp)){
      tmp2 <- tmp[[j]]
      a <- tmp2$formula
      ####add the formula mz
      if(!is.null(a)){
        tmp11<-do.call(rbind,strsplit(as.character(a$Formula),"\\^"))
        tmp12<-sapply(tmp11[,1],function(x){
          DB[,2][which(DB[,1]==x)[1]]
        })
        a<-cbind(a,tmp12)
        form <- paste(paste(a[,1],a[,2],a[,3],sep="_"), collapse=";")
        form<-gsub(" ","",form)
      } else {
        form<-""
      }


      charge <- tmp2$char
      mz <- tmp2$mz
      mz_min <- tmp2$mz_min
      mz_max <- tmp2$mz_max
      RT <- tmp2$rt
      rt_start <- tmp2$rt_start
      rt_end <- tmp2$rt_end
      MS2_lib <- tmp2$consensus
      n_sample <- tmp2$n_sample
      n_iso <- tmp2$n_iso
      total_sample <- length(tmp2$alignment)
      if(!is.null(MS2_lib)){
        #print(i)
        cat(paste("Putative_formula_MS1: ", form, "\n", sep=""), file=fname, append=T)
        cat(paste("Precursor m/z: ", mz, "\n", sep=""), file=fname, append=T)
        cat(paste("Precursor m/z range: ", mz_min,", ", mz_max, "\n", sep=""), file = fname, append=T)
        cat(paste("Charge state: ",charge, "\n", sep=""), file=fname, append=T)
        cat(paste("RT: ", RT, "\n",sep=""),file = fname, append=T)
        cat(paste("RT range: ", rt_start, ", ", rt_end, "\n", sep=""), file=fname, append=T)
        cat(paste("Sample: ", n_sample, "/", total_sample,"\n", sep=""), file=fname, append=T)
        cat(paste("Samples with isotopes: ", n_iso, "\n", sep=""), file=fname, append=T)
        write.table(MS2_lib, fname, sep="\t", row.names=F, col.names=F, append=T, quote=F)
      }
    }
  }
}

