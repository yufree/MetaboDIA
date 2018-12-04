## This function extracts DIA MS2 data using a DDA consensus library
##@param
#### file: Output file name
#### DDA_lib: DDA consensus library object
#### input: A data object containing all the processed information for DIA samples
#### DDA_DIA_compound_map: A list object that maps the formula in DDA library and DIA formula identification
#### com_map: A data object that maps the DIA samples by the identified molecular formula
#############for example, com_map[[1]] <-> formula 1, length(com_map[[1]])=# of samples
#############com_map[[1]][[1]]contain all the index (e.g 300,1001), then input[[300]] and input[[1001]]
#############are all identified as formula 1
#### compound_list_DDA: the formula list in the consensus DDA library
#### f: DIA file names
#### RT_expand: Time span (in seconds) used for alignment (refer to the methods section in the paper)
#### ms1_tol: Mass tolerance for ms1 precursor ion peak (in Dalton, it will be ignored if ms1_ppm is set)
#### ms1_ppm: Mass tolerance for ms1 precursor ion peak (in ppm)
#### ms2_tol: Mass tolerance for ms2 fragment ion peak (in Dalton, it will be ignored if ms2_ppm is set)
#### ms2_ppm: Mass tolerance for ms2 fragment ion peak (in ppm)
#### RT_tol: the RT tolerance when mapping the formula in DDA library with the CIUs in DIA samples
#### cv_n: Minimum number of samples to compute CV from.

DDA_DIA_ms2 <- function(file, DDA_lib, input, DDA_DIA_compound_map, com_map, compound_list_DDA, f, RT_expand, ms1_ppm, ms1_tol, ms2_tol, ms2_ppm, RT_tol, cv_n){

  ####output header
  fields <- c("formula", "charge", "RT", "DDA_lib_mz", "DDA_lib_intensity", "CV", "missing", as.vector(row.names(f)))
  write.table(t(fields), file, sep="\t", row.names=F, col.names=F, quote=F)

  ##store the CVs

  for(y in 1:length(DDA_DIA_compound_map)){
    DIA_id <- DDA_DIA_compound_map[[y]]
    libs <- list()
    if(length(DIA_id) > 0){
      tmp <- DDA_lib[[y]]
      formu <- NULL
      char_rt <- NULL
      extracted <- list()
      for(i in 1:length(tmp)){
        if(!is.null(tmp[[i]]$consensus)){
          libs[[length(libs)+1]] <- tmp[[i]]$consensus
          char_rt <- rbind(char_rt,c(tmp[[i]]$char, tmp[[i]]$rt))
          extracted[[length(extracted)+1]] <- matrix(0,nrow=nrow(tmp[[i]]$consensus), ncol=length(input))
          formu <- c(formu,paste(as.vector(tmp[[i]]$formula[,1]), collapse = ";"))
        }
      }
      if(length(libs)>0){
        DIA_tmp <- com_map[[DIA_id]]

        ####add the monoisotopic peak intensities
        map_rt <- lapply(c(1:length(DIA_tmp)),function(x){
          xx <- DIA_tmp[[x]]
          if(length(xx) > 0){
            sapply(xx,function(yy){
              #yy=xx[1]
              c(yy, input[[x]]$metab[[yy]]$mz_ms1, input[[x]]$metab[[yy]]$charge)
            })
          }
        })
        temp <- do.call(cbind,map_rt)
        charge_l <- unique(temp[3,])
        ref_mz <- sapply(charge_l,function(x){
          c(x,median(temp[2,which(temp[3,]==x)]))
        })
        #####
        for(i5 in 1:length(DIA_tmp)){

          tmp <- input[[i5]]$metab
          mzs <- do.call(rbind, lapply(tmp, function(x) c(x$mz_ms1, x$charge)))
          if(!is.null(ms1_ppm)) ms1_tol <- ref_mz[2,] * ms1_ppm / 1000000

          for(i6 in 1:ncol(ref_mz)){
            tmp_mzs <- abs(mzs[,1] - ref_mz[2,i6])
            id_tmp <- which(mzs[,2] == ref_mz[1,i6] & tmp_mzs<=ms1_tol[i6])
            if(length(id_tmp)>0){
              DIA_tmp[[i5]] <- union(DIA_tmp[[i5]],id_tmp)
            }
          }
        }
        ####end add the mono peak intensities


        map_with_ms2 <- vector("list",length(input))
        for(i4 in 1:length(input)){
          if(!is.null(DIA_tmp[[i4]])){
            DIA_tmp2 <- which(DIA_tmp[[i4]] %in% input[[i4]]$with_ms2)
            if(length(DIA_tmp2) > 0){
              map_with_ms2[[i4]] <- DIA_tmp[[i4]][DIA_tmp2]
            }
          }
        }
        for(i4 in 1:length(input)){
          tp1 <- map_with_ms2[[i4]]
          if(length(tp1) > 0){
            for(j in 1:length(tp1)){
              DIA_char <- input[[i4]]$metab[[tp1[j]]]$charge
              DIA_rt <- input[[i4]]$metab[[tp1[j]]]$RT
              DIA_ms2 <- input[[i4]]$metab[[tp1[j]]]$ms2 ###possible list
              DIA_bp <- input[[i4]]$metab[[tp1[j]]]$basepeak ###possible list
              ###which one should map to
              char_mat <- which(char_rt[,1]==DIA_char)
              if(length(char_mat) > 0){
                RT_diff <- abs(char_rt[char_mat,2] - DIA_rt)
                if(RT_diff[which.min(RT_diff)] < RT_tol){

                  scores <- unlist(lapply(libs,function(x) {
                    ###
                    max(unlist(lapply(DIA_ms2, function(xx){
                      compare_ms2(x,xx,ms2_tol,ms2_ppm)
                    })))

                  }))
                  scores[which(RT_diff>=RT_tol)] <- 0
                  match_to <- which.max(scores)
                  ######found which one matched to and match the MS2
                  if(max(scores)>0){
                    lib_to_match <- libs[[match_to]]
                    tmp_scores <- unlist(lapply(DIA_ms2, function(xxx){
                      compare_ms2(lib_to_match,xxx,ms2_tol,ms2_ppm)
                    }))
                    w <- which.max(tmp_scores)
                    ms2_tmp <- DIA_ms2[[w]]
                    ms2_bp <- DIA_bp[[w]]
                    ms2_tmp[,2] <- ms2_tmp[,2]*ms2_bp
                    ms2 <- ms2_tmp
                    #####match to library and


                    tmp_ints <- sapply(lib_to_match[,1],function(x){
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
                    extracted[[match_to]][,i4] <- apply(cbind(extracted[[match_to]][,i4],tmp_ints),1,max)
                  }# if max score >0
                }#if RT diff<RT_tol
              }# if there is charge match to lib
            }#for each elution
          }#if there is elution with ms2 in the sample
        }#for  each sample
      }#if there is DDA library
    }##if DIA_id>0

    ###output

    if(length(libs)>0){
      for ( i_cv in 1:length(libs)){
        con_lib <- matrix(libs[[i_cv]][which(apply(extracted[[i_cv]],1,sum)>0),],ncol=2)
        ms2_align <- matrix(extracted[[i_cv]][which(apply(extracted[[i_cv]],1,sum)>0),], ncol=length(input))
        #### calculate CV
        cv_score <- apply(ms2_align,1,function(xx){
          tmp_cv <- log2(xx[which(xx>0)])
          if(length(tmp_cv)>=cv_n){
            sd(tmp_cv)/mean(tmp_cv)
          } else {
            return (NA)
          }
        })
        missing <- apply(ms2_align,1,function(xx){
          1-length(which(xx>0))/length(xx)
        })

        out <- cbind(rep(char_rt[i_cv,1], length(cv_score)), rep(char_rt[i_cv,2],length(cv_score)),
                     matrix(con_lib, nrow = length(cv_score)), cv_score, missing,
                     matrix(ms2_align, nrow = length(cv_score)))
        write.table(out, file, append = T, sep="\t", row.names=rep(formu[i_cv],length(cv_score)), col.names=F, quote=F)

      }
    }

    ##output


  }##for DDA_DIA_com_list
}
