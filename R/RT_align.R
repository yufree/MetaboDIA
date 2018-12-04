## Given the same molecular formula, this function RT-aligns the peak features as described in the paper
##@param
####map_charge: a list with lenght of number of samples. Each element ontains the information elution profiles
################from the sample and these elution profiles have the same identified formula and charge
#### RT_expand: Time window (in seconds) for the alignment (refer to the methods section in the paper)


RT_align <- function(map_charge, RT_expand){
  n_sample <- lapply(map_charge, ncol)
  tp=NULL
  for(i in 1:length(n_sample)){
    if(!is.null(n_sample[[i]]))
      tp <- rbind(tp, cbind(i,1:n_sample[[i]]))
  }
  alignment <- list()
  tmp1 <- matrix(do.call(cbind, map_charge)[3:4,], nrow=2)
  tmp1 <- rbind(tmp1, c(1:ncol(tmp1)))

  ###expand RT for alignment
  tmp1[1,] <- tmp1[1,]-RT_expand
  tmp1[2,] <- tmp1[2,]+RT_expand

  while(ncol(tmp1)>0){
    ######reformat the data
    tmp <- rbind(cbind(tmp1[1,], 1, tmp1[3,]),cbind(tmp1[2,], 0, tmp1[3,]))
    tmp <- tmp[order(tmp[,1]), ]

    #####overlapping and segmentation
    seg <- list()
    n <- 0
    seg2elution.i <- vector()
    for(i in 1:(nrow(tmp)-1)){
      if(tmp[i,2]==1){
        n=n+1
        seg2elution.i <- c(seg2elution.i,tmp[i,3])
      } else {
        n=n-1
        seg2elution.i <-  seg2elution.i[!seg2elution.i==tmp[i,3]]
      }
      seg[[i]] <- list(cn=n,seg=seg2elution.i, map_rt=NULL)
    }

    ########output the highest copy number
    copynumber <- unlist(lapply(seg,function(x){x$cn}))
    id <- which(copynumber==max(copynumber))
    rm <- NULL
    for(j in 1:length(id)){
      s <- seg[[id[j]]]$seg
      rm <- union(rm, s)
    }
    ##########check the overlapping between segments only if there are multi max cn
    idl <- as.list(id)
    ####check the overlap with first seg
    while (length(idl)>1){
      ol <- 1
      while(ol==1&length(idl)>1){
        ol <- 0
        for(k in 2:length(idl)){
          if(length(intersect(seg[[idl[[1]]]]$seg, seg[[idl[[k]]]]$seg))>0 ){
            seg[[idl[[1]]]]$seg <- union(seg[[idl[[1]]]]$seg, seg[[idl[[k]]]]$seg)
            idl[[k]] <- NULL
            ol <- 1
            break
          }
        }
      }
      ##output the first seg
      alignment[[length(alignment)+1]] <- seg[[idl[[1]]]]
      idl[[1]] <- NULL
    }
    if(length(idl)==1) alignment[[length(alignment)+1]] <- seg[[idl[[1]]]]




    ######remove the coresopoding illution and restart
    tmp1 <- matrix(tmp1[,-sapply(rm,function(x){which(tmp1[3,]==x)})], nrow=3)

  }
  for(k2 in 1:length(alignment)){
    tmp4 <- matrix(tp[alignment[[k2]]$seg,], ncol=2)
    map_rt <- vector("list", length(map_charge))
    for(l in 1:nrow(tmp4)){
      map_rt[[tmp4[l,1]]] <- matrix(cbind(map_rt[[tmp4[l,1]]], map_charge[[tmp4[l,1]]][ ,tmp4[l,2]]),nrow = 6)
    }
    alignment[[k2]]$map_rt <- map_rt
  }
  alignment

}
