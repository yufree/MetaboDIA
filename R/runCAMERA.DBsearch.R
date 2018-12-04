## This function runs CAMERA and DBsearch
##@param
#### dir: Path to the CAMERA output
#### DB.file: Database file name
#### adduct.file: Adduct file name
#### mode: Positive or negative
#### n_core: Number of cores
#### ppm: Tolerance when doing the DBsearch
#### parameters for CAMERA: prefilter=c(2,50), method="centWave",peakwidth=c(5,60),cor_eic_th = 0.6,scanrange = NULL

runCAMERA.DBsearch<-function(dir=getwd(), DB.file, adduct.file, mode, n_core=4,
                             ppm=10, prefilter=c(2,50), method="centWave", peakwidth=c(5,60), cor_eic_th = 0.6,
                             scanrange = NULL){

  library(CAMERA)
  library(parallel)
  mclapply <- switch( Sys.info()[['sysname']],
                      Windows = {mclapply.hack},
                      Linux   = {mclapply},
                      Darwin  = {mclapply})
  setwd(dir)
  inf<-list.files(pattern = ".mzXML", ignore.case=T)
  outf<-sub(".mzXML",".cam.csv", inf, ignore.case=T)
  outf2<-sub(".mzXML",".cam.metab.csv", inf, ignore.case=T)
  #DB<-read.csv(DB.file)
  DB <- read.delim(DB.file, header=T, sep="\t", as.is=T)
  adduct<-read.csv(adduct.file)
  adduct<-adduct[which(adduct$Use==1),]
  dummy<-mclapply(as.list(c(1:length(inf))), mc.silent=T, mc.cores=n_core, function(x){
    ####run CAMERA
    xs <- xcmsSet(inf[x], method=method, ppm=ppm, peakwidth=peakwidth, prefilter=prefilter, scanrange=scanrange)
    an   <- xsAnnotate(xs, polarity=mode)
    anF  <- groupFWHM(an)
    anI  <- findIsotopes(anF, ppm=ppm)
    anIC <- groupCorr(anI, cor_eic_th=cor_eic_th)
    peaklist <- getPeaklist(anIC)
    write.csv(peaklist, file=outf[x], quote=F)

    ###conver the CAMERA output and run DBsearch
    Q<-convert.camera.output(inf[x], mode)
    Q_mass<-lapply(Q$Q1,function(y) y$mz )
    Q_charge<-lapply(Q$Q1,function(y) y$Charge )
    form_list<-lapply(1:length(Q_mass),function(zz) {
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
            paste(DB$formula[id[[zzz]]], adduct$Ion.name[ind[zzz]], sep="^")
          }
        }))
        best_form<-formu[which.min(unlist(mz_diff))]
        all_form<-unlist(lapply(1:length(id), function(zzz) {
          if(length(id[[zzz]])>0){
            paste(unique(DB$formula[id[[zzz]]]), adduct$Ion.name[ind[zzz]], sep="^")
          }
        }))
        return(c(paste(all_form,collapse = ";"), best_form))
      } else {
        return(c(NA,NA))
      }
    })
    out<-do.call(rbind, form_list)
    tmp<-data.frame(index.Q=c(1:length(form_list)), hit.DB=out[,1], best.hit=out[,2])
    write.csv(tmp,outf2[x], row.names=F, quote=F)
  })
}
