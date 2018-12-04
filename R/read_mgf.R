## This function reads the MS2 information form mgf file (For DIA_umpire output)
##@param
####fname: file name
####min_spectrum: minimum number of MS2 peaks (otherwise the MS2 epctrum will be discarded)

read_mgf <- function(fname, min_spectrum){
  mgf <- scan(fname,what="", sep="$", quiet=T)
  #locate information
  begin <- grep("BEGIN", mgf)
  end <- grep("END", mgf)
  scan <- grep("scan", mgf)
  RT <- grep("RTINSECONDS", mgf)
  pre <- grep("PEPMASS", mgf)
  charge <- grep("CHARGE", mgf)
  if(length(charge)>0){
    b <- as.vector(sapply(sub("CHARGE=","", mgf[charge]), switch,
                        "1+" = 1,
                        "2+" = 2,
                        "3+" = 3,
                        "4+" = 4,
                        "5+" = 5,
                        "6+" = 6,
                        "1-" = -1,
                        "2-" = -2,
                        "3-" = -3,
                        "4-" = -4,
                        "5-" = -5,
                        "6-" = -6))
  } else {
    b <- NULL
  }
  ###number of lines before the spectrum data
  n_fields <- grep("^[[:digit:]]", mgf[1:1000])[1]-begin[1]
  dat <- cbind(begin+n_fields,end-1)
  RT <- as.numeric(sub("RTINSECONDS=","",mgf[RT]))
  pre <- as.numeric(sub("PEPMASS=","",mgf[pre]))
  a <- list()
  for(i in 1:nrow(dat)){
    a[[i]] <- dat[i,]
  }
  d <- lapply(a,function(x){
    if(x[2]-x[1]>0){
      tmp <- do.call(rbind, lapply(as.list(mgf[x[1]:x[2]]), function(y){
        as.numeric(unlist(strsplit(y," ")))
      }))

      mz <- tmp[,1]
      count <- tmp[,2]
      ms2 <- cbind(mz,count)
      ms2
    } else {
      ms2 <- NA
      ms2
    }
  })

  out <- list()
  for(i in 1:length(RT)){
    if(!is.na(d[[i]])[1]){
      if(nrow(d[[i]])>min_spectrum){
        out[[length(out)+1]] <- list(RT=RT[i], ms1=pre[i], charge=b[i], ms2=d[[i]])
      }
    }
  }
  out
}
