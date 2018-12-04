## This function is MS1 only workflow: equivalent to Part 1 in the DDA_DIA_workflow
## without processing the MS2 information here

MS1_extraction <- function(file_dir, n_core=4, ms1_output="MS1_extraction.txt", RT_expand=6, cv_n=5,
                         ms1_ppm=10, ms1_tol=0.001, metmatch_FDR=NA, best_hit=F, intensity_type="area"){
  library(parallel)
  mclapply  <-  switch( Sys.info()[['sysname']],
                      Windows = {mclapply.hack},
                      Linux   = {mclapply},
                      Darwin  = {mclapply})
  setwd(file_dir)
  fl <- list.files(pattern=".cam.csv")
  sample_name <- gsub(".cam.csv", "", fl)
  dbsearch_name <- gsub(".cam.csv", ".cam.metab.csv", fl)
  iso_file <- gsub(".cam.csv", ".cam.iso.csv", fl)
  mono_file <- gsub(".cam.csv", ".cam.mono.csv", fl)
  cutoff <- gsub(".cam.csv", ".cam.cutoffs.txt", fl)
  file_check <- which(file.exists(dbsearch_name)=="FALSE")
  if(length(file_check)>0){
    stop(paste(paste(dbsearch_name[file_check], collapse = ", "), "do/does not exist. Please check file_matching.txt ",sep=" "))
  }
  file_output <- data.frame(sample=sample_name, MS1_elution=fl, dbsearch=dbsearch_name)
  write.table(file_output, "file_matching.txt", sep="\t", row.names=F, quote=F)

  f <- read.table("file_matching.txt", header=T, sep="\t", row.names=1)
  input_DDA_files <- list()

  for(i in 1:nrow(f)){
    input_DDA_files[[i]] <- list(um_out=as.character(iso_file[i]), dbsearch=as.character(dbsearch_name[i]),
                               cutoffs=cutoff[i], monopeak=mono_file[i])
  }
  input_DDA <- mclapply(input_DDA_files, mc.silent=T, mc.cores = n_core, function(x){
    readmetab_DDA_2(x[[1]], x[[2]], x[[3]], metmatch_FDR, best_hit, intensity_type)
  })

  clist_DDA <- mclapply(input_DDA, mc.silent=T, mc.cores=n_core, function(x){
    tmp <- lapply(x$metab,function(y) unlist(strsplit(y$DBcompound,";")) )
  })
  ###############build
  if(best_hit){
    compound_list_DDA <- unique(unlist(lapply(clist_DDA,unlist  )))
    compound_list_DDA <- compound_list_DDA[which(!is.na(compound_list_DDA))]
    compound_map_DDA <- list()

    for(i in 1:length(input_DDA_files)){
      a <- clist_DDA[[i]]
      compound_map_DDA[[i]] <- mclapply(as.list(compound_list_DDA), mc.silent=T, mc.cores=n_core, function(x){
        which(a==x)
      })
    }
    com_map_DDA <- list()
    for(i in 1:length(compound_list_DDA)){
      tmp <- list()
      for(j in 1:length(input_DDA_files)){
        tmp[[j]] <- compound_map_DDA[[j]][[i]]
      }
      com_map_DDA[[i]] <- tmp
    }

  } else {
    clist_DDA_2 <- clist_DDA
    tmp_list <- do.call(c, clist_DDA_2)
    tl <- unlist(lapply(tmp_list, length))
    tmp_id1 <- which(tl>1)
    tmp_id2 <- which(tl==1)
    compound_list_DDA <- list()
    compound_list_2 <- NULL
    #aa <- NULL
    for(i in tmp_id1){

      tmp_10 <- tmp_list[[i]]
      if(length(intersect(compound_list_2, tmp_10))==0){
        compound_list_2 <- c(compound_list_2, tmp_10)
        compound_list_DDA[[length(compound_list_DDA)+1]] <- tmp_10
      } else {
        compound_list_2 <- union(compound_list_2, tmp_10)
        tmp_id10 <- unlist(lapply(compound_list_DDA, function(xx){
          length(intersect(xx, tmp_10))
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
      if(!is.element(tmp_10, compound_list_2)){
        compound_list_2 <- c(compound_list_2, tmp_10)
        compound_list_DDA[[length(compound_list_DDA)+1]] <- tmp_10
      }
    }

    compound_map_DDA <- mclapply(as.list(c(1:length(input_DDA_files))), mc.silent=T, mc.cores=n_core, function(i){
      a <- clist_DDA_2[[i]]
      abc <- lapply(as.list(compound_list_DDA), function(x){
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

  }

  mono_map <- mclapply(input_DDA_files, mc.cores=n_core, mc.silent=T, function(x){
    read_mono_2(mono_peak_file=x[[4]], ms1_tol,ms1_ppm, intensity_type)
  })

  ### add to "com_map_DDA"
  for(x in 1:length(mono_map)){
    input_DDA[[x]]$with_ms2 <- c(input_DDA[[x]]$with_ms2, (mono_map[[x]]$with_ms2+length(input_DDA[[x]]$metab)))
    input_DDA[[x]]$metab <- c(input_DDA[[x]]$metab, mono_map[[x]]$metab)
  }

  #####process each compount map and build library


  #####build library

  DDA_lib <- mclapply(com_map_DDA, mc.silent=T , mc.cores=n_core, function(x){
    com_consensus_lib <- list()


    ####alignment regardless of ms2
    #check RT isobars and charge
    map_rt_charge <- list()
    for(i in 1:length(input_DDA_files)){
      map_rt_charge[[i]] <- sapply(x[[i]],function(y){
        #c(y, input_DDA[[i]]$metab[[y]]$RT, input_DDA[[i]]$metab[[y]]$RT_range, input_DDA[[i]]$metab[[y]]$charge)
        c(y, input_DDA[[i]]$metab[[y]]$RT, input_DDA[[i]]$metab[[y]]$RT_range,
          input_DDA[[i]]$metab[[y]]$charge, input_DDA[[i]]$metab[[y]]$mz_ms1)
      })
    }
    charge_list <- lapply(map_rt_charge,function(z){
      if(length(z)>0) z[5,]
    })


    ######seperate different charges and RT isobars
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
      #######
      for(i3 in 1:length(alignment_rt)){

        ##########
        map_rt <- alignment_rt[[i3]]$map_rt
        temp <- do.call(cbind,map_rt)
        rt_mean <- mean(temp[2,])
        ###find a reference peak
        rt_lower <- min(temp[3,]) - RT_expand
        rt_upper <- max(temp[4,]) + RT_expand
        ref_sample <- min(which(lapply(map_rt,is.null)==FALSE))
        check_id <- which(lapply(map_rt,is.null)==TRUE)
        ref_mz <- input_DDA[[ref_sample]]$metab[[ map_rt[[ref_sample]][1,1] ]]$mz_ms1
        for(i5 in check_id){
          tmp <- input_DDA[[i5]]$metab
          mzs <- do.call(rbind,lapply(tmp,function(x) c(x$mz_ms1, x$RT, x$RT_range, x$charge)))
          if(!is.null(ms1_ppm)) ms1_tol <- ref_mz*ms1_ppm/1000000
          mzs2 <- abs(mzs[,1]-ref_mz)
          id_tmp <- which( mzs2<=ms1_tol & mzs[,3]<= rt_upper & mzs[,4] >= rt_lower)
          if(length(id_tmp)>0){
            map_rt[[i5]] <- matrix(t(cbind(id_tmp, matrix(mzs[id_tmp,c(2,3,4,5,1)], ncol=5))), nrow=6)
          }

        }
        com_consensus_lib[[length(com_consensus_lib)+1]] <- list(char=charge, consensus=NULL,
                                                                 cn=alignment_rt[[i3]]$cn, rt=rt_mean, alignment=map_rt)
      }
    }
    com_consensus_lib
  })

  out_ms1(ms1_output, DDA_lib, input_DDA, compound_list_DDA, f, cv_n)
}


