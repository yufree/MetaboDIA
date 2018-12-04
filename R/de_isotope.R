##### This function removes the isotopes in consensus library from the DIA_DIA_workfolw.
##@params:
#### consen_lib: The consensus library built from DIA_DIA_workflow
#### PROTON_MASS: PROTON_MASS = 1.007276
#### de_ppm: Mass tolerance for detecting isotope peaks

de_isotope<-function(consen_lib,PROTON_MASS,de_ppm){
  if(is.null(consen_lib)){
    return(consen_lib)
  } else {
    if(nrow(consen_lib)<2) {
      return (consen_lib)
    } else{
      ####check isotope peaks from the largest m/z
      iso_check<-sapply(2:nrow(consen_lib),function(x){
        if(1000000*abs((consen_lib[x,1]-(consen_lib[x-1,1]+PROTON_MASS)))/consen_lib[x,1]<de_ppm & consen_lib[x,2]/consen_lib[x-1,2]<1) {
          c1<-1
        } else {
          c1<-0
        }
        if(1000000*abs((consen_lib[x,1]-(consen_lib[x-1,1]+0.5*PROTON_MASS)))/consen_lib[x,1]<de_ppm & consen_lib[x,2]/consen_lib[x-1,2]<1) {
          c2<-1
        } else {
          c2<-0
        }
        max(c1,c2)
      })
      iso_check<-c(0,iso_check)

      return(matrix(consen_lib[which(iso_check==0),], ncol=2))
    }
  }
}
