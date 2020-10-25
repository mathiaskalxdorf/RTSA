dataNormalizationofcond2 = function(cfg_info,data_list,norm_in_log2=T) ###Normalization of condition2 relative to condition1
{
  num_temps <- length(which(grepl("^Temp[0-9]$",colnames(cfg_info$cfg_sheet))))
  Correctiontablecond2 <- as.data.frame(matrix(ncol=(3*cfg_info$num_samples)+1,nrow=num_temps))
  colnames(Correctiontablecond2) <- c("Temp",paste("Densitymax_",cfg_info$cfg_sheet$MSExperimentIDs[1:cfg_info$num_samples],sep=""),paste("Correctionfactor_",cfg_info$cfg_sheet$MSExperimentIDs[1:cfg_info$num_samples],sep=""),paste("NumProts_",cfg_info$cfg_sheet$MSExperimentIDs[1:cfg_info$num_samples],sep=""))
  Correctiontablecond2$Temp <- 1:num_temps
  
  pdf(paste("Normalization of ",cfg_info$condition_name_2,".pdf",sep=""))
  for(g in 1:cfg_info$num_gradientparts) #### go through data of specific gradient parts
  {
    dataindices <- which(cfg_info$cfg_sheet$Gradient == g)
    for(t in 1:num_temps)
    {
      ####get cond1 columns for respective temp in respective experiment
      column <- NULL
      columnnotnorm <- NULL
      columnreffc <- NULL
      colindex <- NULL
      colindexnotnorm <- NULL
      colindexreffc <- NULL
      for(i in dataindices)
      {
        column <- append(column,paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[t,i],"_norm",sep=""))
        columnnotnorm <- append(columnnotnorm,paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[t,i],sep=""))
        columnreffc <- append(columnreffc,paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t,i],"_norm",sep=""))
        colindex <- append(colindex,which(colnames(data_list[[i]]) == paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[t,i],"_norm",sep="")))
        colindexnotnorm <- append(colindexnotnorm,which(colnames(data_list[[i]]) == paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[t,i],sep="")))
        colindexreffc <- append(colindexreffc,which(colnames(data_list[[i]]) == paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t,i],"_norm",sep="")))
      }
      ###extract max densities for adjusted following density plots
      maxdens <- NA
      densymax <- 0###for density plots
      for(i in dataindices)
      {
        ###filter for good quantified proteins
        temp <- data_list[[i]]
        temp <- subset(temp,qusm >= cfg_info$qusm_filter_norm & qupm >= cfg_info$qupm_filter_norm)
        temp <- temp[which(temp[,colindex[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
        temp <- temp[which(temp[,colindexreffc[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
        ####calculate relative fc to respective condition1
        if(norm_in_log2 == 1)
        {
          foldchange <- log2(as.numeric(temp[,colindex[which(dataindices == i)]]) / as.numeric(temp[,colindexreffc[which(dataindices == i)]]))
        }else
        {
          foldchange <- as.numeric(temp[,colindex[which(dataindices == i)]]) / as.numeric(temp[,colindexreffc[which(dataindices == i)]])
        }
        
        Correctiontablecond2[t,1+i] <- maxDensity(foldchange)
        Correctiontablecond2[t,1+(2*cfg_info$num_samples)+i] <- nrow(temp)
        if(is.na(maxdens))
        {
          maxdens <- Correctiontablecond2[t,1+i]
        }else
        {
          maxdens <- append(maxdens,Correctiontablecond2[t,1+i])
        }
        if(length(foldchange[!is.na(foldchange)]) > 0)
        {
          density <- density(foldchange[!is.na(foldchange)])
          if(max(density$y) > densymax){densymax <- max(density$y)}
        }
      }
      if(densymax != 0)
      {
        ####plot density distribution before normalization
        cols <- c("black","red","blue","darkgreen","darkgrey","chocolate","burlywood4","cadetblue","chocolate4","chartreuse4","coral4","deeppink3") ### 12 colors for up to 12 replicates
        ind <- 0
        for(i in dataindices)
        {
          temp <- data_list[[i]]
          temp <- subset(temp,qusm >= cfg_info$qusm_filter_norm & qupm >= cfg_info$qupm_filter_norm)
          temp <- temp[which(temp[,colindex[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
          temp <- temp[which(temp[,colindexreffc[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
          ind <- ind + 1
          if(norm_in_log2 == 1)
          {
            foldchange <- log2(as.numeric(temp[,colindex[which(dataindices == i)]]) / as.numeric(temp[,colindexreffc[which(dataindices == i)]]))
          }else
          {
            foldchange <- as.numeric(temp[,colindex[which(dataindices == i)]]) / as.numeric(temp[,colindexreffc[which(dataindices == i)]])
          }
          
          density <- density(foldchange[!is.na(foldchange)])
          
          if(ind == 1)
          {
            if(norm_in_log2 == 1)
            {
              plot(density,xlim=c(-2,2),ylim=c(0,densymax),col=cols[ind],main=paste("Rel.fcs. (log2) in ",cfg_info$condition_name_2," at ",cfg_info$temp_gradient[cfg_info$temp_order[t,g]]," C",sep=""))
            }else
            {
              plot(density,xlim=c(0,2),ylim=c(0,densymax),col=cols[ind],main=paste("Rel.fcs. in ",cfg_info$condition_name_2," at ",cfg_info$temp_gradient[cfg_info$temp_order[t,g]]," C",sep=""))
            }
          }else
          {
            lines(density,col=cols[ind])
          }
          segments(density$x[which(density$y==max(density$y))], max(density$y), density$x[which(density$y==max(density$y))], par("usr")[3],lty=3, col = cols[ind])
        }
        legend("topright",col=cols[1:ind],legend = paste("Replicate ",1:length(dataindices),sep=""),lty=1)
        
        ####normalize
        for(i in dataindices)
        {
          if(norm_in_log2 == 1)
          {
            corrfactor <- 1 / 2^Correctiontablecond2[t,1+i]#1/maxdensity
          }else
          {
            corrfactor <- 1 / Correctiontablecond2[t,1+i]#1/maxdensity
          }
          Correctiontablecond2[t,1+cfg_info$num_samples+i] <- corrfactor
          normcolname <- column[which(dataindices == i)] #rel.fc
          data_list[[i]][,normcolname] <- data_list[[i]][,colindex[which(dataindices == i)]] * corrfactor #rel.fc.
        }
        
        ####plot density distribution after normalization
        cols <- c("black","red","blue","darkgreen","darkgrey","chocolate","burlywood4","cadetblue","chocolate4","chartreuse4","coral4","deeppink3") ### 12 colors for up to 12 replicates
        ind <- 0
        for(i in dataindices)
        {
          temp <- data_list[[i]]
          temp <- subset(temp,qusm >= cfg_info$qusm_filter_norm & qupm >= cfg_info$qupm_filter_norm)
          temp <- temp[which(temp[,colindexnotnorm[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
          temp <- temp[which(temp[,colindexreffc[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
          ind <- ind + 1
          #normcolname2 <- column2 #sumionarea
          if(norm_in_log2 == 1)
          {
            foldchange <- log2(as.numeric(temp[,colindex[which(dataindices == i)]]) / as.numeric(temp[,colindexreffc[which(dataindices == i)]]))
          }else
          {
            foldchange <- as.numeric(temp[,colindex[which(dataindices == i)]]) / as.numeric(temp[,colindexreffc[which(dataindices == i)]])
          }
          density <- density(foldchange[!is.na(foldchange)])
          
          if(ind == 1)
          {
            if(norm_in_log2 == 1)
            {
              plot(density,xlim=c(-2,2),ylim=c(0,densymax),col=cols[ind],main=paste("Normalized rel.fcs. (log2) in ",cfg_info$condition_name_2," at ",cfg_info$temp_gradient[cfg_info$temp_order[t,g]]," C",sep=""))
            }else
            {
              plot(density,xlim=c(0,2),ylim=c(0,densymax),col=cols[ind],main=paste("Normalized rel.fcs. in ",cfg_info$condition_name_2," at ",cfg_info$temp_gradient[cfg_info$temp_order[t,g]]," C",sep=""))
            }
          }else
          {
            lines(density,col=cols[ind])
          }
          segments(density$x[which(density$y==max(density$y))], max(density$y), density$x[which(density$y==max(density$y))], par("usr")[3],lty=3, col = cols[ind])
        }
        legend("topright",col=cols[1:ind],legend = paste("Replicate ",1:length(dataindices),sep=""),lty=1)
      }else
      {
        Correctiontablecond2[t,1+cfg_info$num_samples+i] <- NA
        normcolname <- column[which(dataindices == i)] #rel.fc
        data_list[[i]][,normcolname] <- data_list[[i]][,colindex[which(dataindices == i)]]
      }
      
    }
  }
  dev.off()
  
  output <- list()
  output[["data_list_norm"]] <- data_list
  output[["cond2_normalization_parameter_sheet"]] <- Correctiontablecond2
  return(output)
  
}