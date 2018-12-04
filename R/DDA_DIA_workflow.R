## This is the main wrapper function for the DDA_DIA_workflow.


##@params

## DIA_dir: The path to the directory containing DIA data
## DDA_dia: The path to the directory containing DDA data
## n_core: # of cores used for mclapply

###MS2 consensus library critiria
## min_sample=0.2: At least 20% of samples has MS2 data to build MS2 consensus library
## consensus_filter=0.05: Average relative intensity cutoff
## ms2_rep=0.5: At least 50% reproducibility across the MS2 spectra
##              for a specific ms2 peak in the consensus library
## ms2_IQR=80: The inter quartile range of the mz differences between
##             the consensus peak and the peaks that align to that peak
##             is 80 ppm,  which is 2*ms2_ppm
## ms2_ppm=40: The mass tolerance in ppm when aligning consensus MS2 peaks
## ms2_tol=0.01: The mass tolerance in Dalton when aligning consensus MS2 peaks (ignored when ms2_ppm is set)
## The library is reported in "consensus_library.txt"

###RT alignment
## RT_expand=6: we expand the XIC elution time for +/-6 seconds when XICs are aligned across the samples

###Linking MS2 spectrum and MS1 profile
## ms1_ppm=10: The mass tolerance in ppm when link the MS1 peak m/z and the precursor ion m/z in MS2 data
## ms2_tol=0.01 the mass tolerance in Dalton it is ignored when ms1_ppm is set

###Read and clean MS2 data
## min_noise=20: Noise level in MS2 data
## min_int=0.01: The minimum relative intensity of MS2 data
## min_basepeak=100: The minimum intensity of the basepeak
## min_spectrum=2: The minimum number of samples with ms2 spectra


###Output option
## ms1_output: File name for MS1 intensity alignment
## ms2_output: File name for MS2 extraction for DIA_extraction
## file: File name for MS2 extraction for DDA_DIA_extraction
## plot_lib_var=TRUE plot the consensus library quality the m/z difference in ppm  file name is "consensus_library_quality.pdf"

###MetMatch options (not used in this version)
## metmatch_FDR=NA
## best_hit=FALSE

###Additional formula identification option
## DB.file: Database file name
## adduct.file: adduct file name
## identification.ppm: The mass tolerance in ppm when doing the additional formula identification

DDA_DIA_workflow <- function(DIA_dir, DDA_dir, mode, n_core=4, min_sample=0.5, consensus_filter=0.05,
                           ms2_rep=0.5, ms2_IQR=80, RT_expand=6, r_score=0.5,
                           ms1_ppm=10, ms1_tol=0.001, ms2_ppm=40, ms2_tol=0.01,
                           min_noise=20, min_int=0.01, min_basepeak =100, min_spectrum=1,
                           RT_tol=60,cv_n=5, ms2_mz_tol=0.5,identification.ppm=10,
                           ###added features: Use metmatch or not, run with decoy, FDR cutoff, use best.hit or all hits
                           metmatch_FDR=NA, best_hit=F, intensity_type="area", DB.file, adduct.file){

  library(parallel)
  mclapply  <-  switch( Sys.info()[['sysname']],
                      Windows = {mclapply.hack},
                      Linux   = {mclapply},
                      Darwin  = {mclapply})
  PROTON_MASS = 1.007276 # Mass of Proton
  file="DDA_DIA_extraction.txt"

  ########## PART 1: DDA library building
  print("Processing DDA data......")

  setwd(DDA_dir)

  ####read all the information
  fl1 <- list.files(pattern=".mzxml", ignore.case = T)
  fl <- gsub(".mzxml", ".cam.csv", fl1, ignore.case = T)
  sample_name <- gsub(".cam.csv","",fl)
  dbsearch_name <- gsub(".cam.csv", ".cam.metab.csv", fl)
  file_check <- which(file.exists(dbsearch_name) == "FALSE")
  iso_file <- gsub(".cam.csv", ".cam.iso.csv", fl)
  mono_file <- gsub(".cam.csv", ".cam.mono.csv", fl)
  cutoff <- gsub(".cam.csv", ".cam.cutoffs.txt", fl)
  if(length(file_check)>0){
    stop(paste(paste(dbsearch_name[file_check],collapse = ", "),"do/does not exist. Please check file_matching.txt ",sep=" "))
  }
  ms2_name <- fl1
  file_output <- data.frame(sample=sample_name,MS1_elution=fl,dbsearch=dbsearch_name)
  write.table(file_output,"file_matching.txt",sep="\t",row.names = F, quote=F)

  #####################
  f_DDA <- read.table("file_matching.txt",header=T,sep="\t",row.names = 1)
  input_DDA_files <- list()

  for(i in 1:length(fl)){
    input_DDA_files[[i]] <- list(um_out=as.character(iso_file[i]),dbsearch=as.character(dbsearch_name[i]),
                               ms2=as.character(ms2_name[i]),cutoffs=cutoff[i],monopeak=mono_file[i])
  }

  ####process the peaks with isotope peaks
  input_DDA <- mclapply(input_DDA_files,mc.cores = n_core,mc.silent=T,function(x){
    readmetab_DDA(x[[1]],x[[2]],x[[3]],x[[4]],
                  ms1_tol, ms1_ppm, min_spectrum, min_basepeak,
                  min_int, min_noise, mode, ms2_mz_tol,
                  metmatch_FDR, best_hit, intensity_type)
  })

  ####process the peaks without isotope peaks
  mono_map <- mclapply(input_DDA_files,mc.cores =n_core,mc.silent=T,function(x){
   read_mono(mgf_file = x[[3]], mono_peak_file = x[[5]],
             ms1_tol,ms1_ppm, min_spectrum, min_basepeak,
             min_int, min_noise, mode, ms2_mz_tol, intensity_type)
  })


  ##addition formula identification steps for the monoisotpe peaks
  #DB <- read.csv(DB.file)
  DB <- read.delim(DB.file, header=T, sep="\t", as.is=T)
  adduct <- read.csv(adduct.file)
  adduct <- adduct[which(adduct$Use==1),]
  mono_2 <- mclapply(as.list(c(1:length(input_DDA_files))), mc.cores = n_core, function(x){
    metab <- NULL
    with_ms2 <- NULL
    tmp <- mono_map[[x]]$metab
    for(i in 1:length(tmp)){
      tmp2 <- mono_map[[x]]$metab[[i]]
      ####### note this filter applies when build library and remove when MS1 extraction
      if(length(tmp2$ms2)>0){
        z <- tmp2$mz_ms1
        c.i <- tmp2$charge
        ####search for different adducts
        ind <- which(adduct$Charge==c.i)
        mz_converted <- (z-adduct$adduct.mass[ind])*adduct$Mult[ind]

        id <- lapply(mz_converted,function(zzz) which(abs(DB$exactMass-zzz)/zzz <= identification.ppm/1000000))

        if(length(unlist(id))>0){
          mz_diff <- lapply(mz_converted, function(zzz) abs(DB$exactMass[which(abs(DB$exactMass-zzz)/zzz <= identification.ppm/1000000)]-zzz))
          formu <- unlist(lapply(1:length(id), function(zzz) {
            if(length(id[[zzz]])>0){
              paste(DB$formula[id[[zzz]]], adduct$Ion.name[ind[zzz]], sep="^")
            }
          }))
          best_form <- formu[which.min(unlist(mz_diff))]
          all_form <- unlist(lapply(1:length(id),function(zzz) {
            if(length(id[[zzz]])>0){
              paste(unique(DB$formula[id[[zzz]]]), adduct$Ion.name[ind[zzz]], sep="^")
            }
          }))
        } else {
          all_form <- NA
          best_form <- NA
        }
        if(!is.na(best_form)){
          tmp3 <-  mono_map[[x]]$metab[[i]]
          tmp3$formula.all <-  paste(all_form,collapse = ";")
          tmp3$DBcompound <-  best_form
          metab[[length(metab)+1]] <- tmp3
        }


      }
    }
    with_ms2 <- c(1:length(metab))
    list(metab=metab,with_ms2=with_ms2)
  })
  mono_map <-  mono_2
  rm(mono_2)


  ###some data structure transformation
  input_2 <- mclapply(as.list(c(1:length(input_DDA_files))), mc.cores = n_core, function(x){
    metab <- NULL
    with_ms2 <- NULL
    tmp <- input_DDA[[x]]$metab
    for(i in 1:length(tmp)){
      tmp2 <- input_DDA[[x]]$metab[[i]]
      if(length(tmp2$ms2)>0){
        tmp3 <-  input_DDA[[x]]$metab[[i]]
        metab[[length(metab)+1]] <- tmp3
      }
    }
    with_ms2 <- c(1:length(metab))
    list(metab=metab,with_ms2=with_ms2)
  })

  input_DDA <- input_2
  rm(input_2)

  ### add mono_map to "com_map_DDA"
  for(x in 1:length(mono_map)){
    input_DDA[[x]]$metab <- c(input_DDA[[x]]$metab,mono_map[[x]]$metab)
    input_DDA[[x]]$with_ms2 <- c(1:length(input_DDA[[x]]$metab))
  }

  ### identify the formula list and build a map for each formula
  clist_DDA <- mclapply(input_DDA,mc.silent = T,mc.cores=n_core,function(x){
    tmp <- lapply(x$metab,function(y) unlist(strsplit(y$DBcompound,";")) )
  })

  clist_DDA_2 <-  clist_DDA
  tmp_list <- do.call(c, clist_DDA_2)
  tl <- unlist(lapply(tmp_list,length))
  tmp_id1 <- which(tl>1)
  tmp_id2 <- which(tl==1)
  compound_list_DDA <- list()
  compound_list_2 <- NULL
  for(i in tmp_id1){

    tmp_10 <- tmp_list[[i]]
    if(length(intersect(compound_list_2,tmp_10))==0){
      compound_list_2 <- c(compound_list_2,tmp_10)
      compound_list_DDA[[length(compound_list_DDA)+1]] <- tmp_10
    } else {
      compound_list_2 <- union(compound_list_2,tmp_10)
      tmp_id10 <- unlist(lapply(compound_list_DDA,function(xx){
        length(intersect(xx,tmp_10))
      }))
      id_11 <- which(tmp_id10>0)
      tmp_11 <- tmp_10
      for(j in length(id_11):1){
        tmp_11 <- union(tmp_11,compound_list_DDA[[id_11[j]]])
        compound_list_DDA[[id_11[j]]] <- NULL
      }
      compound_list_DDA[[length(compound_list_DDA)+1]] <- tmp_11

    }
  }

  for(i in tmp_id2){
    tmp_10 <- tmp_list[[i]]
    if(!is.element(tmp_10,compound_list_2)){
      compound_list_2 <- c(compound_list_2,tmp_10)
      compound_list_DDA[[length(compound_list_DDA)+1]] <- tmp_10
    }
  }

  compound_map_DDA <- mclapply(as.list(c(1:length(input_DDA_files))), mc.silent = T, mc.cores=n_core, function(i){
    a <- clist_DDA_2[[i]]
    abc <- lapply(as.list(compound_list_DDA),function(x){
      tmp_check <- lapply(a,function(y) length(intersect(x,y))>0)
      which(tmp_check==T)
    })
    abc
  })
  com_map_DDA <- list()
  for(i in 1:length(compound_list_DDA)){
    tmp <- list()
    for(j in 1:length(input_DDA_files)){
      tmp[[j]] <- compound_map_DDA[[j]][[i]]
    }
    com_map_DDA[[i]] <- tmp
  }

  #####identify CIUs and build DDA based consensus library

  print("Building DDA consensus library......")
  DDA_lib <- mclapply(as.list(c(1:length(com_map_DDA))), mc.silent=T, mc.cores=n_core, function(ii){
    #check if have ms2 info
    x <- com_map_DDA[[ii]]
    com_consensus_lib <- list()
    ####alignment regardless of ms2


    #check RT isobars and charge
    map_rt_charge <- list()
    for(i in 1:length(input_DDA_files)){
      map_rt_charge[[i]] <- sapply(x[[i]],function(y){
        c(y,input_DDA[[i]]$metab[[y]]$RT, input_DDA[[i]]$metab[[y]]$RT_range,
          input_DDA[[i]]$metab[[y]]$charge, input_DDA[[i]]$metab[[y]]$mz_ms1)
      })
    }
    charge_list <- lapply(map_rt_charge,function(z){
      if(length(z)>0) z[5,]
    })

    ######seperate different charges only
    charge_list <- unique(unlist(charge_list))
    for(i in 1:length(charge_list)){
      charge <- charge_list[i]
      map_charge <- lapply(map_rt_charge,function(x){
        if(is.matrix(x)){
          if(length(which(x[5,]==charge))>0) as.matrix(x[,which(x[5,]==charge)])
        }
      })
      ######seperate RT isobars

      alignment_rt <- RT_align( map_charge,RT_expand)
      #collecting the ms2 spectrum
      for(i3 in 1:length(alignment_rt)){

        ##########
        map_rt <- alignment_rt[[i3]]$map_rt
        mz_median <- NA
        mz_min <- NA
        mz_max <- NA
        rt_mean <- NA
        rt_start <- NA
        rt_end <- NA
        n_iso <- 0

        map_with_ms2 <- map_rt
        formulae <- NULL
        MS1_quant <- NULL
        mean_elution_time <- NULL
        #####check if have ms2 spectrum
        if(sum(unlist(lapply(map_with_ms2,length)))>0){



          ####start here library build recursion function
          recursion_check <- 1
          tmp_input_DDA <- input_DDA

          while(recursion_check==1){
            n_sample <- 0
            n_spec <- 0
            n_iso <- 0
            ms2_spectrum <- list()

            for(i2 in 1:length(input_DDA_files)){
              if(!is.null(map_with_ms2[[i2]])&length(map_with_ms2[[i2]])>0){
                n_sample <- n_sample+1
                id <- map_with_ms2[[i2]][1,]
                iso_increment <- 0
                for(j2 in 1:length(id)){
                  tmp <- tmp_input_DDA[[i2]]$metab[[ id[j2] ]]$ms2
                  if(input_DDA[[i2]]$metab[[ id[j2] ]]$N_iso>1) iso_increment <- 1
                  for (k2 in 1:length(tmp)){
                    ms2_spectrum[[length(ms2_spectrum)+1]] <- tmp[[k2]]
                  }
                }
                n_iso <- n_iso+iso_increment
              }
            }

            #########build consensus library
            if(n_sample>=length(input_DDA_files)*min_sample){
              temp <- lapply(map_with_ms2,function(x) {
                if(length(x)>0) x
              })
              temp <- do.call(cbind,temp)
              mz_median <- median(temp[6,])
              mz_min <- min(temp[6,])
              mz_max <- max(temp[6,])
              rt_mean <- median(temp[2,])
              rt_start <- min(temp[3,])
              rt_end <- max(temp[4,])
              temp_lib <- consensus(ms2_spectrum, ms2_tol, ms2_ppm, consensus_filter, ms2_rep,ms2_IQR)
              consen_lib <- temp_lib$lib
              mz_diff <- temp_lib$mz_diff
              n_spec <- length(ms2_spectrum)
              if(!is.null(consen_lib)){ ## consensus lib built, check the min qc_score
                qc_score <- lapply(ms2_spectrum,function(y){
                  compare_ms2(consen_lib, y, ms2_tol, ms2_ppm)
                })
                ms2_spectra <- ms2_spectrum
                if(min(unlist(qc_score))< r_score){ ###remove the low score MS2
                  ####check every MS2 spectrum


                  for(i_check in 1:length(input_DDA_files)){
                    if(!is.null(map_with_ms2[[i_check]])&length(map_with_ms2[[i_check]])>0){##each sample
                      id_check <- map_with_ms2[[i_check]][1,]
                      elution_to_keep <- rep(1,length(id_check)) ##if multi elution which ones to keep
                      for(j_check in 1:length(id_check)){
                        tmp_check <- tmp_input_DDA[[i_check]]$metab[[ id_check[j_check] ]]$ms2
                        check_score <- unlist(lapply(tmp_check,function(xx){
                          compare_ms2(consen_lib,xx,ms2_tol,ms2_ppm)
                        }))
                        high_score_id <- which(check_score >= r_score)
                        if(length(high_score_id)>0){
                          tmp_list <- NULL
                          for(i2_check in 1:length(high_score_id)){
                            tmp_list[[length(tmp_list)+1]] <- tmp_input_DDA[[i_check]]$metab[[ id_check[j_check] ]]$ms2[[i2_check]]
                          }
                          tmp_input_DDA[[i_check]]$metab[[ id_check[j_check] ]]$ms2 <- tmp_list
                          rm(tmp_list)
                        } else {
                          elution_to_keep[j_check] <- 0
                        }
                      } ##end for j_check
                      if(length(which(elution_to_keep>0))==0) {
                        map_with_ms2[[i_check]] <- list()
                      } else {
                        map_with_ms2[[i_check]] <- matrix(map_with_ms2[[i_check]][,which(elution_to_keep>0)], nrow=6)
                      }

                    }
                  }
                } else { ## all pass: stop the checking
                  recursion_check <- 0
                }

              } else { ###no consensus lib built: exit while loop
                qc_score <- NULL
                n_spec <- 0
                ms2_spectra <- NULL
                recursion_check <- 0
              }

            } else {
              consen_lib <- NULL
              mz_diff <- NULL
              qc_score <- NULL
              ms2_spectra <- NULL
              n_spec <- 0
              recursion_check <- 0
            }
          }    ####end while(check==1)
          rm(tmp_input_DDA)


        } else {
          n_sample <- 0
          n_spec <- 0
          n_iso <- 0
          consen_lib <- NULL
          mz_diff <- NULL
          qc_score <- NULL
          ms2_spectra <- NULL
        } ##end if else there is ms2 spectrum
        ####
        ###record  MS1 quantification and formula
        if(!is.null(consen_lib)){
          formulae <- table(unlist(lapply(1:length(map_with_ms2),function(xx){
            x <- map_with_ms2[[xx]]
            if(!is.null(x)&length(x)>0){
              unique(unlist(sapply(x[1,],function(y){
                if(!is.na(input_DDA[[xx]]$metab[[y]]$DBcompound)){
                  strsplit( input_DDA[[xx]]$metab[[y]]$formula.all,";")[[1]]
                }

              })))
            }
          })))
          if(length(formulae)>0){
            formulae <- data.frame(Formula=names(formulae), Count=as.vector(formulae))
            formulae <- formulae[order(formulae$Count),]
          } else {
            formulae <- data.frame(Formula=compound_list_DDA[[ii]], Count=NA)
            formulae <- formulae[order(formulae$Count, decreasing = T), ]
          }

          #####  MS1 quantification
          MS1_quant <- unlist(lapply(1:length(map_with_ms2), function(xx){
            x <- map_with_ms2[[xx]]
            if(!is.null(x)&length(x)>0){
              tmp_int <- sapply(x[1,], function(y){
                tmp_ms1_int <- input_DDA[[xx]]$metab[[y]]$intensity_ms1
                tmp_check <- input_DDA[[xx]]$metab[[y]]$ms2
                check_score <- max(unlist(lapply(tmp_check, function(z){
                  compare_ms2(consen_lib,z, ms2_tol, ms2_ppm)
                })))
                c(tmp_ms1_int,check_score)
              })
              return(max(tmp_int[1,which(tmp_int[2,] == max(tmp_int[2,]))]))
            } else {
              return (NA)
            }


          }))
          ###MS1 mean elution time
          mean_elution_time <- unlist(lapply(1:length(map_with_ms2),function(xx){
            x <- map_with_ms2[[xx]]
            if(!is.null(x)&length(x)>0){
              tmp_int <- sapply(x[1,],function(y){
                tmp_ms1_int <- input_DDA[[xx]]$metab[[y]]$intensity_ms1
                tmp_rt <- input_DDA[[xx]]$metab[[y]]$RT_range
                tmp_check <- input_DDA[[xx]]$metab[[y]]$ms2
                check_score <- max(unlist(lapply(tmp_check,function(z){
                  compare_ms2(consen_lib,z,ms2_tol,ms2_ppm)
                })))
                c(tmp_ms1_int,check_score,tmp_rt[2]-tmp_rt[1])
              })
              tt <- matrix(tmp_int[ ,which(tmp_int[2,]==max(tmp_int[2,]))], nrow=3)
              return(max(tt[3,which(tt[1, ]==max(tt[1, ]))]))
            } else {
              return (NA)
            }
          }))
          mean_elution_time <- mean(mean_elution_time, na.rm = T)

        }


        com_consensus_lib[[length(com_consensus_lib)+1]] <- list(formula=formulae,MS1_quant=MS1_quant,mean_elution_time=mean_elution_time,char=charge,consensus=consen_lib,rt=rt_mean,rt_start=rt_start,rt_end=rt_end,
                                                               mz=mz_median, mz_min=mz_min, mz_max=mz_max,mz_var=mz_diff,QC=qc_score,ms2=ms2_spectra,
                                                               n_sample=n_sample,n_iso=n_iso,n_spec=n_spec,cn=alignment_rt[[i3]]$cn,alignment=map_with_ms2)

      }
    }
    com_consensus_lib
  })


  ###output MS1 intensitis based on the consensus library
  bb <- mclapply(DDA_lib,mc.cores = n_core,function(x){

    aa <- lapply(x,function(y){
      if(!is.null(y$consensus)){
        a <- y$formula
        char <- y$char
        rt <- y$rt
        b <- y$MS1_quant
        ref_mz<-y$mz
        cc <- y$mean_elution_time
        a <- paste(paste(a[,1], a[,2], sep="_"),collapse=";")
        miss <- length(which(is.na(b))) / length(b)
        tmp_cv <- log2(b[which(b>0)])
        if(length(tmp_cv) >= cv_n){
          cv <- sd(tmp_cv) / mean(tmp_cv)
        } else {
          cv <- NA
        }
        c(a,ref_mz,char,rt,cv,miss,cc,b)
      } else {
        return(NULL)
      }
    })
    do.call(rbind,aa)
  })
  ms1out <- do.call(rbind,bb)
  fields <- c("formula","ref_mz","charge","RT","CV","missing","mean_elution_time", as.vector(row.names(f_DDA)))
  write.table(t(fields), "Lib_MS1_quant.txt", sep="\t",row.names=F, col.names=F, quote=F)

  write.table(ms1out, "Lib_MS1_quant.txt", sep="\t", row.names=F, col.names=F, append = T, quote=F)


  ####output consensus library
  out_library(DDA_lib,file="DDA_library.txt",DB)



  ##############Part 2: process DIA samples the process is same as DDA part
  print("Processing DIA data......")
  setwd(DIA_dir)
  fl <- list.files(pattern=".cam.csv")
  sample_name <- gsub(".cam.csv","",fl)
  dbsearch_name <- gsub(".cam.csv",".cam.metab.csv",fl)
  iso_file <- gsub(".cam.csv",".cam.iso.csv",fl)
  mono_file <- gsub(".cam.csv",".cam.mono.csv",fl)
  cutoff <- gsub(".cam.csv",".cam.cutoffs.txt",fl)
  file_check <- which(file.exists(dbsearch_name)=="FALSE")
  if(length(file_check)>0){
    stop(paste(paste(dbsearch_name[file_check],collapse = ", "),"do/does not exist. Please check file_matching.txt ",sep=" "))
  }
  ms2_name1 <- gsub(".cam.csv","_Q1.mgf",fl)
  file_check <- which(file.exists(ms2_name1)=="FALSE")
  if(length(file_check)>0){
    stop(paste(paste(ms2_name1[file_check],collapse = ", "),"do/does not exist. Please check file_matching.txt ",sep=" "))
  }
  ms2_name2 <- gsub(".cam.csv","_Q2.mgf",fl)
  file_check <- which(file.exists(ms2_name2)=="FALSE")
  if(length(file_check)>0){
    stop(paste(paste(ms2_name2[file_check],collapse = ", "),"do/does not exist. Please check file_matching.txt ",sep=" "))
  }
  file_output <- data.frame(sample=sample_name, MS1_elution=fl, dbsearch=dbsearch_name, MS2_Q1=ms2_name1, MS2_Q2=ms2_name2)
  write.table(file_output,"file_matching.txt",sep="\t",row.names = F, quote=F)


  ##################################
  f <- read.table("file_matching.txt",header=T,sep="\t",row.names=1)
  input_files <- list()

  for(i in 1:length(iso_file)){
    input_files[[i]] <- list(um_out=as.character(iso_file[i]),dbsearch=as.character(dbsearch_name[i]),
                           Q1=as.character(ms2_name1[i]),Q2=as.character(ms2_name2[i]), cutoffs=cutoff[i], monopeak=mono_file[i])
  }

  input <- mclapply(input_files,mc.silent = T,mc.cores=n_core, function(x){
    readmetab_DIA(x[[1]],x[[2]],x[[3]],x[[4]],x[[5]],
                  ms1_tol, ms1_ppm, min_spectrum, min_basepeak,
                  min_int, min_noise, mode, ms2_mz_tol,
                  metmatch_FDR, best_hit, intensity_type)
  })
  mono_map <- mclapply(input_files, mc.cores=n_core, mc.silent=T, function(x){
    read_mono(mgf_file = c(x[[3]],x[[4]]), mono_peak_file = x[[6]],
              ms1_tol, ms1_ppm, min_spectrum, min_basepeak,
              min_int, min_noise, mode, ms2_mz_tol, intensity_type)
  })
  #DB <- read.csv(DB.file)
  DB <- read.delim(DB.file, header=T, sep="\t", as.is=T)
  mono_2 <- mclapply(as.list(c(1:length(input_files))), mc.cores = n_core, function(x){
    metab <- NULL
    with_ms2 <- NULL
    tmp <- mono_map[[x]]$metab
    for(i in 1:length(tmp)){
      tmp2 <- mono_map[[x]]$metab[[i]]
      ####### note this filter applies when build library and remove when MS1 extraction
      if(length(tmp2$ms2)>0){
        z <- tmp2$mz_ms1
        c.i <- tmp2$charge
        ####search for different adducts
        ind <- which(adduct$Charge == c.i)
        mz_converted <- (z-adduct$adduct.mass[ind]) * adduct$Mult[ind]

        id <- lapply(mz_converted,function(zzz) which(abs(DB$exactMass-zzz)/zzz <= identification.ppm / 1000000))

        if(length(unlist(id))>0){
          mz_diff <- lapply(mz_converted,function(zzz) abs(DB$exactMass[which(abs(DB$exactMass-zzz)/zzz <= identification.ppm / 1000000)]-zzz))
          formu <- unlist(lapply(1:length(id),function(zzz) {
            if(length(id[[zzz]])>0){
              paste(DB$formula[id[[zzz]]], adduct$Ion.name[ind[zzz]], sep="^")
            }
          }))
          best_form <- formu[which.min(unlist(mz_diff))]
          all_form <- unlist(lapply(1:length(id),function(zzz) {
            if(length(id[[zzz]])>0){
              paste(unique(DB$formula[id[[zzz]]]), adduct$Ion.name[ind[zzz]], sep="^")
            }
          }))
        } else {
          all_form <- NA
          best_form <- NA
        }
        if(!is.na(best_form)){
          tmp3 <-  mono_map[[x]]$metab[[i]]
          tmp3$formula.all <-  paste(all_form, collapse = ";")
          tmp3$DBcompound <-  best_form
          metab[[length(metab)+1]] <- tmp3
        }
      }
    }
    with_ms2 <- c(1:length(metab))
    list(metab=metab, with_ms2=with_ms2)
  })
  mono_map <-  mono_2
  rm(mono_2)
  input_2 <- mclapply(as.list(c(1:length(input_files))), mc.cores = n_core, function(x){
    metab <- NULL
    with_ms2 <- NULL
    tmp <- input[[x]]$metab
    for(i in 1:length(tmp)){
      tmp2 <- input[[x]]$metab[[i]]
      ####### note this filter applies when build library and remove when MS1 extraction
      if(length(tmp2$ms2)>0){
        tmp3 <-  input[[x]]$metab[[i]]
        metab[[length(metab)+1]] <- tmp3
      }
    }
    with_ms2 <- c(1:length(metab))
    list(metab=metab, with_ms2=with_ms2)
  })

  input <- input_2
  rm(input_2)
  ### add to "com_map_DDA"
  for(x in 1:length(mono_map)){
    input[[x]]$metab <- c(input[[x]]$metab, mono_map[[x]]$metab)
    input[[x]]$with_ms2 <- c(1:length(input[[x]]$metab))
  }

  ####unique compound list
  clist_DIA <- mclapply(input, mc.silent = T, mc.cores=n_core, function(x){
    tmp <- lapply(x$metab, function(y) unlist(strsplit(y$DBcompound,";")) )
  })

  ###############build
  clist_DIA_2 <-  clist_DIA
  tmp_list <- do.call(c, clist_DIA_2)
  tl <- unlist(lapply(tmp_list, length))
  tmp_id1 <- which(tl>1)
  tmp_id2 <- which(tl==1)
  compound_list_DIA <- list()
  compound_list_2 <- NULL
  for(i in tmp_id1){

    tmp_10 <- tmp_list[[i]]
    if(length(intersect(compound_list_2,tmp_10))==0){
      compound_list_2 <- c(compound_list_2,tmp_10)
      compound_list_DIA[[length(compound_list_DIA)+1]] <- tmp_10
    } else {
      compound_list_2 <- union(compound_list_2,tmp_10)
      tmp_id10 <- unlist(lapply(compound_list_DIA,function(xx){
        length(intersect(xx,tmp_10))
      }))
      id_11 <- which(tmp_id10>0)
      id_11 <- id_11[order(id_11)]
      tmp_11 <- tmp_10
      for(j in length(id_11):1){
        tmp_11 <- union(tmp_11,compound_list_DIA[[id_11[j]]])
        compound_list_DIA[[id_11[j]]] <- NULL
      }
      compound_list_DIA[[length(compound_list_DIA)+1]] <- tmp_11

    }
  }

  for(i in tmp_id2){
    tmp_10 <- tmp_list[[i]]
    if(!is.element(tmp_10,compound_list_2)){
      compound_list_2 <- c(compound_list_2,tmp_10)
      compound_list_DIA[[length(compound_list_DIA)+1]] <- tmp_10
    }
  }

  ######Part 3: using DDA based library
  ##map the formulae between DDA library and DIA data
  print("Extracting DIA MS/MS intensities using DDA consensus library......")
  id_to_rm <- NULL
  form_to_add <- list()
  for(i in 1:length(compound_list_DDA)){
    tm1 <-  compound_list_DDA[[i]]
    if(length(tm1)>1){ ##merge possible compound_list_DIA
      tmp_check <- lapply(compound_list_DIA,function(y) length(intersect(tm1,y))>0)
      tm_id <- which(tmp_check==T)
      if(length(tm_id)>1){
        id_to_rm <- c(id_to_rm,tm_id)
        tm_compound <- unique(as.vector(unlist(sapply(tm_id,function(z) compound_list_DIA[[z]]))))
        form_to_add[[length(form_to_add)+1]] <- tm_compound
      }

    }
  }
  form_to_add_2 <- list()
  while(length(form_to_add)>0){
    ck <- lapply(form_to_add, function(xx) length(intersect(xx,form_to_add[[1]])) )
    ck_id <- which(ck>0)
    if(length(ck_id>1)){
      ck_id <- ck_id[order(ck_id,decreasing = T)]
      form_to_add_2[[length(form_to_add_2)+1]] <- unique(unlist(lapply(as.list(ck_id),function(x) form_to_add[[x]])))
      for(i in ck_id) form_to_add[[i]] <- NULL
    } else {
      fom_to_add_2[[length(fom_to_add_2)+1]] <- form_to_add[[1]]
      fom_to_add[[1]] <- NULL
    }

  }
  if(length(id_to_rm)>0){
    id_to_rm <- id_to_rm[order(id_to_rm,decreasing = T)]
    for(i in id_to_rm) compound_list_DIA[[i]] <- NULL
  }
  if(length(form_to_add_2)>0){
    for(i in 1:length(form_to_add_2)) compound_list_DIA[[length(compound_list_DIA)+1]] <-  form_to_add_2[[i]]
  }


  ###############build
  compound_map <- mclapply(as.list(c(1:length(input_files))), mc.silent = T, mc.cores=n_core, function(i){
    a <- clist_DIA_2[[i]]
    abc <- lapply(as.list(compound_list_DIA),function(x){
      tmp_check <- lapply(a,function(y) length(intersect(x,y))>0)
      which(tmp_check==T)
    })
    abc
  })

  ##transpose the list for the next step
  com_map <- list()
  for(i in 1:length(compound_list_DIA)){
    tmp <- list()
    for(j in 1:length(input_files)){
      tmp[[j]] <- compound_map[[j]][[i]]
    }
    com_map[[i]] <- tmp
  }

  DDA_DIA_compound_map <- lapply(compound_list_DDA,function(x){
    tmp_check <- lapply(compound_list_DIA,function(y) length(intersect(x,y))>0)
    which(tmp_check==T)
  })

  ######output DIA MS2 information
  DDA_DIA_ms2(file, DDA_lib, input, DDA_DIA_compound_map, com_map, compound_list_DDA,
              f, RT_expand, ms1_ppm, ms1_tol, ms2_tol, ms2_ppm, RT_tol, cv_n)

}

