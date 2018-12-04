## This function reads the camera output and converts to the data structure 
## required for the ensuing process.
## It also seperates the output into two files: peak features with and without isotopes
## @param
#### fh: file name of camera output
#### mode: positive or negative (ionization mode)

convert.camera.output<-function (fh,mode)
{
  d = read.csv(fh)
  e = d[grepl("M", d$isotopes) == FALSE, ]
  e = e[order(e$mz), ]

  d = d[grep("M", d$isotopes), ]
  d = d[order(d$isotopes), ]
  isotopes = as.character(d$isotopes)
  id = unlist(lapply(isotopes, function(x) strsplit(x, "\\[")[[1]][2]))
  id = as.numeric(sub("\\]", "", id))
  u.id = unique(id)
  d = cbind(d, id)
  Q1 = list()
  isoto <- NULL
  for (i in u.id) {
    pattern = paste("^", i, "$", sep = "")
    tm = d[grepl(pattern, d$id), ]
    M = tm$mz
    if(mode=="positive"){
      charge = round(1/(M[2]-M[1]), 0)
    } 
    else {
      charge = -round(1/(M[2]-M[1]), 0)
    }

    isoto <- rbind(isoto,cbind(tm, charge=charge))
    Q1[[i]] = list(Charge=charge, RT = tm$rt[1], mz = M[1], Nisotope = length(M) )
  }
  write.csv(isoto,sub(".csv",".iso.csv",fh), quote=FALSE)
  write.csv(e,sub(".csv",".mono.csv",fh), quote=FALSE)

  list(Q1=Q1)
}
