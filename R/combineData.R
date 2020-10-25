addCols = function(appenddataframe,originaldataframe,colsname) ###used for reordering columns for combined data table
{
  for(col in colsname)
  {
    columns <- which(grepl(col,colnames(originaldataframe)))
    appenddataframe <- cbind(appenddataframe,originaldataframe[,columns])
  }
  return(appenddataframe)
}

combineData = function(summary_data_list_normalized,cfg_info,fordata = "both") ###combine data for data analysis
{
  
  if(fordata == "both")
  {
    calcsetup <- c(1,2)
  }else
  {
    if(fordata == "with")
    {
      calcsetup <- c(1)
    }
    if(fordata == "without")
    {
      calcsetup <- c(2)
    }
  }
  summary_data_combined <- list()
  
  for(ndata in calcsetup)
  {
    combined <- summary_data_list_normalized[[ndata]][[1]]
    
    for(i in 2:cfg_info$num_samples)
    {
      combined <- full_join(combined,summary_data_list_normalized[[ndata]][[i]],by=c("gene_name"="gene_name"))
    }
    
    ###Reorder columns
    combinedreorder <- as.data.frame(combined[,1])###Protein
    colnames(combinedreorder) <- "Protein"
    
    ####general data
    combinedreorder <- addCols(appenddataframe = combinedreorder,originaldataframe = combined,colsname = c("representative","proteinscore","totalpsm","qusm\\.","qupm\\.","ms1intensity","ms1seqs","ms1maxminusmin"))
    
    ####relfcs for each temp
    for(t in cfg_info$temp_gradient)
    {
      combinedreorder <- addCols(appenddataframe = combinedreorder,originaldataframe = combined,colsname = c(paste(t,"C_cond1",sep=""),paste(t,"C_cond2",sep="")))
    }
    ####log2ratios for each temp
    for(t in cfg_info$temp_gradient)
    {
      combinedreorder <- addCols(appenddataframe = combinedreorder,originaldataframe = combined,colsname = paste("log2Ratio",t,sep=""))
    }
    #####calculate meanlog2ratios
    for(t in cfg_info$temp_gradient)
    {
      cols <- which(grepl(paste("log2Ratio",t,sep=""),colnames(combinedreorder)))
      mean <- as.data.frame(rowMeans(combinedreorder[,cols],na.rm=T))
      mean[is.na(mean)] <- NA
      colnames(mean) <- paste("meanlog2Ratio",t,"C",sep="")
      combinedreorder <- cbind(combinedreorder,mean)
    }
    #####calculate stdevs based on log2ratios
    stdevs <- as.data.frame(matrix(ncol=length(cfg_info$temp_gradient),nrow=nrow(combinedreorder)))
    for(t in cfg_info$temp_gradient)
    {
      tindex <- which(cfg_info$temp_gradient == t)
      cols <- which(grepl(paste("^log2Ratio",t,"C",sep=""),colnames(combinedreorder)))
      stdevs[,tindex] <- apply(combinedreorder[,cols],1,function(x) sd(x,na.rm=T))
      colnames(stdevs)[tindex] <- paste("StDev_",t,"C",sep="")
    }
    stdev <- apply(stdevs,1,function(x) mean(x,na.rm=T))###average SD at all temps 
    stdevs <- cbind(stdevs,stdev)
    colnames(stdevs)[ncol(stdevs)] <- "meanStDev"
    stdevs[is.na(stdevs)] <- NA
    combinedreorder <- cbind(combinedreorder,stdevs)
    if(ndata == 1){summary_data_combined[["with_abundance_effect"]] <- combinedreorder}
    if(ndata == 2){summary_data_combined[["without_abundance_effect"]] <- combinedreorder}
  }
  return(summary_data_combined)
}
