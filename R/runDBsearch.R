
## This function runs DBsearch for molecular formula assignment
##@param
#### dir: Path to the CAMERA output
#### DB.file: Database file name
#### adduct.file: Adduct file name
#### mode: Positive or negative
#### n_core: Number of cores
#### ppm: Mass tolerance for DBsearch


runDBsearch<-function(dir=getwd(), DB.file, adduct.file, mode, n_core=4, ppm=10){
  library(parallel)
  mclapply <- switch( Sys.info()[['sysname']],
                      Windows = {mclapply.hack},
                      Linux   = {mclapply},
                      Darwin  = {mclapply})
  setwd(dir)
  inf<-list.files(pattern = ".cam.csv")
  outf2<-sub(".cam.csv",".cam.metab.csv",inf)
  #DB<-read.csv(DB.file)
  DB <- read.delim(DB.file, header=T, sep="\t", as.is=T)
  adduct<-read.csv(adduct.file)
  adduct<-adduct[which(adduct$Use==1),]
  #adduct
  dummy<-mclapply(as.list(c(1:length(inf))), mc.silent=T, mc.cores=n_core, function(x){
    Q<-convert.camera.output(inf[x], mode)
    Q_mass<-lapply(Q$Q1, function(y) y$mz )
    Q_charge<-lapply(Q$Q1, function(y) y$Charge )
    form_list<-lapply(1:length(Q_mass), function(zz) {
      z<-Q_mass[[zz]]
      c.i<-Q_charge[[zz]]
      ####search for different adducts
      ind<-which(adduct$Charge==c.i)
      mz_converted<-(z-adduct$adduct.mass[ind])*adduct$Mult[ind]

      id<-lapply(mz_converted, function(zzz) which(abs(DB$exactMass-zzz)/zzz<=ppm/1000000))

      if(length(unlist(id))>0){
        mz_diff<-lapply(mz_converted, function(zzz) abs(DB$exactMass[which(abs(DB$exactMass-zzz)/zzz<=ppm/1000000)]-zzz))
        formu<-unlist(lapply(1:length(id), function(zzz) {
          if(length(id[[zzz]])>0){
            paste(DB$formula[id[[zzz]]],adduct$Ion.name[ind[zzz]],sep="^")
          }
        }))
        best_form<-formu[which.min(unlist(mz_diff))]
        all_form<-unlist(lapply(1:length(id), function(zzz) {
          if(length(id[[zzz]])>0){
            paste(unique(DB$formula[id[[zzz]]]),adduct$Ion.name[ind[zzz]],sep="^")
          }
        }))
        return(c(paste(all_form, collapse = ";"), best_form))
      } else {
        return(c(NA,NA))
      }
    })
    out<-do.call(rbind, form_list)
    tmp<-data.frame(index.Q=c(1:length(form_list)), hit.DB=out[,1], best.hit=out[,2])
    write.csv(tmp,outf2[x], row.names=F, quote=F)
  })
}
