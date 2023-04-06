calculateRatios = function(summary_data_list_normalized,cfg_info,fordata = "both")###Calculate ratios between condition2 and condition1
{
  num_temps <- length(which(grepl("^Temp[0-9]$",colnames(cfg_info$cfg_sheet))))
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
  
  for(ndata in calcsetup)
  {
    data_list <- summary_data_list_normalized[[ndata]]
    ###Calculate ratios
    for(i in 1:cfg_info$num_samples)
    {
      for(t in 1:num_temps)
      {
        column <- paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[t,i],"_norm",sep="")
        columnref <- paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t,i],"_norm",sep="")
        ###get col indices
        colindex <- which(colnames(data_list[[i]]) == column)#rel.fc condition2
        colindexref <- which(colnames(data_list[[i]]) == columnref)#rel.fc condition1
        ratiocolname <- paste("Ratio",cfg_info$temp_gradient[cfg_info$temp_order[t,cfg_info$dataset_to_gradientpart$Gradientpart[i]]],"C",sep="")
        ####replace relfc with value = 0 by 0.01 and then calculate ratios
        for(p in 1:nrow(data_list[[i]]))
        {
          if(any(!is.na(data_list[[i]][p,colindexref])))
          {
            if(any(data_list[[i]][p,colindexref] == 0))
            {
              data_list[[i]][p,colindexref] = 0.01
            }
          }
          if(any(!is.na(data_list[[i]][p,colindexref])) & any(!is.na(data_list[[i]][p,colindex])))
          {
            if(any(data_list[[i]][p,colindex] >= cfg_info$relfc_cutoff_plot) | any(data_list[[i]][p,colindexref] >= cfg_info$relfc_cutoff_plot))
            {
              data_list[[i]][p,ratiocolname] <- data_list[[i]][p,colindex]/data_list[[i]][p,colindexref]
            }else
            {
              data_list[[i]][p,ratiocolname] <- NA
            }
          }else
          {
            data_list[[i]][p,ratiocolname] <- NA
          }
        }
        log2ratiocolname <- paste("log2",ratiocolname,sep="")
        data_list[[i]][,log2ratiocolname] <- log2(data_list[[i]][,ratiocolname])
      }
      ####if a protein is present more than once (different IPI IDs) select the protein with the higher protein score
      Data <- data_list[[i]]
      allproteins <- Data$gene_name[!duplicated(Data[,"gene_name"])]
      if(length(allproteins) < nrow(Data)) ###if there are duplicated names
      {
        for(prots in Data$gene_name[which(duplicated(Data[,"gene_name"]))])
        {
          inde <- which(Data$gene_name==prots)
          
          if(length(inde) > 0)
          {
            sub <- subset(Data,gene_name==prots)
            ####keep the entry with the best protein score and remove all other from the analysis
            keepentry <- which(sub$proteinscore==max(sub$proteinscore))[1]
            removeindes <- inde[-keepentry]
            if(length(removeindes) > 0)
            {
              Data <- Data[-removeindes,]
            }
          }
        }
      }
      data_list[[i]] <- Data
    }
    
    ####rename colnames so that they can be later distinguished
    for(i in 1:cfg_info$num_samples)
    {
      colnames(data_list[[i]])[-1] <- paste(colnames(data_list[[i]])[-1],".",cfg_info$dataset_to_gradientpart$Replicate[i],".",cfg_info$dataset_to_gradientpart$Gradientpart[i],sep="")
    }
    
    #####rename rel_fcs based on temperature
    for(g in unique(cfg_info$dataset_to_gradientpart$Gradientpart))
    {
      dataListindices <- which(cfg_info$dataset_to_gradientpart$Gradientpart == g)
      for(i in dataListindices)
      {
        for(t in 1:num_temps)
        {
          colnames(data_list[[i]]) <- gsub(paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t,i],"_norm",sep=""),paste(cfg_info$temp_gradient[cfg_info$temp_order[t,g]],"C_cond1",sep=""),colnames(data_list[[i]]))
          colnames(data_list[[i]]) <- gsub(paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[t,i],"_norm",sep=""),paste(cfg_info$temp_gradient[cfg_info$temp_order[t,g]],"C_cond2",sep=""),colnames(data_list[[i]]))
        }
      }
    }
    
    #####rename sumionareas based on temperature
    for(g in unique(cfg_info$dataset_to_gradientpart$Gradientpart))
    {
      dataListindices <- which(cfg_info$dataset_to_gradientpart$Gradientpart == g)
      for(i in dataListindices)
      {
        for(t in 1:num_temps)
        {
          colnames(data_list[[i]]) <- gsub(paste("sumionarea_protein_",cfg_info$tmt_channel_per_exp_condition1[t,i],"_norm",sep=""),paste("sumionarea_",cfg_info$temp_gradient[cfg_info$temp_order[t,g]],"C__cond1",sep=""),colnames(data_list[[i]]))
          colnames(data_list[[i]]) <- gsub(paste("sumionarea_protein_",cfg_info$tmt_channel_per_exp_condition2[t,i],"_norm",sep=""),paste("sumionarea_",cfg_info$temp_gradient[cfg_info$temp_order[t,g]],"C__cond2",sep=""),colnames(data_list[[i]]))
        }
      }
    }
    
    summary_data_list_normalized[[ndata]] <- data_list
  }
  
  return(summary_data_list_normalized)
}