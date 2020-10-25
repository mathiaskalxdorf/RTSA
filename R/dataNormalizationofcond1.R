dataNormalizationofcond1 = function(cfg_info,data_list,norm_in_log2=T) ###Normalization of all condition1 to each other at each individual temperature based on rel.fcs
{
  ####based on proteins with qusm >= filter criteria, qupm >= filtercriteria and relfcs > 0.4
  #####Save all correction factors for condition1 in one table
  num_temps <- length(which(grepl("^Temp[0-9]$",colnames(cfg_info$cfg_sheet))))
  Correctiontable <- as.data.frame(matrix(ncol=(3*cfg_info$num_samples)+1,nrow=num_temps))
  colnames(Correctiontable) <- c("Temp",paste("Densitymax_",cfg_info$cfg_sheet$MSExperimentIDs[1:cfg_info$num_samples],sep=""),paste("Correctionfactor_",cfg_info$cfg_sheet$MSExperimentIDs[1:cfg_info$num_samples],sep=""),paste("NumProts_",cfg_info$cfg_sheet$MSExperimentIDs[1:cfg_info$num_samples],sep=""))
  Correctiontable$Temp <- 1:nrow(Correctiontable)
  Correctiontable[1,2:(cfg_info$num_samples+1)] <- NA
  Correctiontable[1,(2*cfg_info$num_samples+2):((3*cfg_info$num_samples)+1)] <- NA
  ####Add columns for rel_fc_protein_norm
  for(i in 1:cfg_info$num_samples)
  {
    relfc_col_names <- colnames(data_list[[i]])[which(grepl("^rel_fc_protein_",colnames(data_list[[i]])))]
    
    temp <- data.frame(temp=data_list[[i]][,relfc_col_names])
    colnames(temp) <- paste(relfc_col_names,"_norm",sep="")
    data_list[[i]] <- cbind(data_list[[i]],temp)
  }
  
  ####Normalize condition1 at each temperature between replicates
  pdf(paste("Normalization of ",cfg_info$condition_name_1,".pdf",sep=""))
  for(g in 1:cfg_info$num_gradientparts) #### go through data of specific gradient parts
  {
    dataindices <- which(cfg_info$cfg_sheet$Gradient == g)
    for(t in 1:(num_temps-1)) #### x temperatures are normalized per TMT experiment, reference temp is not normalized
    {
      ####get cond1 columns for respective temp in respective experiment
      column <- NULL
      columnnotnorm <- NULL
      colindex <- NULL
      colindexnotnorm <- NULL
      for(i in dataindices)
      {
        column <- append(column,paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t+1,i],"_norm",sep=""))
        columnnotnorm <- append(columnnotnorm,paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t+1,i],sep=""))
        colindex <- append(colindex,which(colnames(data_list[[i]]) == paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t+1,i],"_norm",sep="")))
        colindexnotnorm <- append(colindexnotnorm,which(colnames(data_list[[i]]) == paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition1[t+1,i],sep="")))
      }
      ###extract max densities
      maxdens <- NA
      densymax <- 0###for density plots
      for(i in dataindices)
      {
        temp <- data_list[[i]]
        temp <- subset(temp,qusm >= cfg_info$qusm_filter_norm & qupm >= cfg_info$qupm_filter_norm)
        temp <- temp[which(temp[,colindex[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
        
        if(norm_in_log2 == 1){Correctiontable[t+1,1+i] <- maxDensity(log2(temp[,colindex[which(dataindices == i)]]))}else{Correctiontable[t+1,1+i] <- maxDensity(temp[,colindex[which(dataindices == i)]])}
        Correctiontable[t+1,1+(2*cfg_info$num_samples)+i] <- nrow(temp)
        if(is.na(maxdens))
        {
          maxdens <- Correctiontable[t+1,1+i]
        }else
        {
          maxdens <- append(maxdens,Correctiontable[t+1,1+i])
        }
        if(nrow(temp)>0)
        {
          if(norm_in_log2 == 1){density <- density(log2(temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]]))}else{density <- density(temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]])}
          if(max(density$y) > densymax){densymax <- max(density$y)}
        }
        
      }
      #####check if there is an outlier in the selected maxdensity and if yes check if there is another possible maxdensity closer to the others
      combinations <- combn(maxdens, length(maxdens)-1)
      combinations2 <- matrix(ncol=ncol(combinations),nrow=nrow(combinations))
      SDpercombination <- colSds(combinations)
      SDall <- sd(maxdens)
      for(m in maxdens) ###convert values to indices
      {
        combinations2[which(combinations == m)] <- which(maxdens == m)
      }
      ###is there a combination removing one maxdens which has a much better SD
      ####SD of combination has to have at least 2x less SD
      indx <- which(2*SDpercombination <= SDall)
      dontnormalize <- F
      if(!is.na(median(SDpercombination)) & median(SDpercombination) > 0.25 & length(indx) == 0)
      {
        if(length(SDpercombination) > 3 & min(SDpercombination) < median(SDpercombination)-0.03)###if there are many replicates then just use the combination with best SD
        {
          indx <- which(SDpercombination == min(SDpercombination))
        }else ###dont normalize at all
        {
          dontnormalize <- T
        }
      }
      if(length(indx) > 0)
      {
        ###which datapoint has to be kicked out to get a much better sd
        range <- 1:length(maxdens)
        outlierindx <- which(range %not in% combinations2[,indx])
        ####now have a look again in that outlier and detect if more than one maxima would be possible
        i <- dataindices[outlierindx]
        temp <- data_list[[i]]
        temp <- subset(temp,qusm >= cfg_info$qusm_filter_norm & qupm >= cfg_info$qupm_filter_norm)
        temp <- temp[which(temp[,colindex[outlierindx]] >= cfg_info$relfc_cutoff_norm),]
        
        if(norm_in_log2 == 1)
        {
          density <- density(log2(temp[!is.na(temp[,colindex[outlierindx]]),colindex[outlierindx]]))
          maxima <- allmaxDensity(vals = log2(temp[!is.na(temp[,colindex[outlierindx]]),colindex[outlierindx]]))
        }else
        {
          density <- density(temp[!is.na(temp[,colindex[outlierindx]]),colindex[outlierindx]])
          maxima <- allmaxDensity(vals = temp[!is.na(temp[,colindex[outlierindx]]),colindex[outlierindx]])
        }
        
        
        if(length(maxima) > 1)
        {
          sds <- matrix(ncol=1,nrow=length(maxima))
          for(m in 1:length(maxima))####go through all possible maxima and check if SD gets reduced
          {
            sds[m] <- sd(append(maxdens[-outlierindx],maxima[m]))
          }
          indximproved <- which(sds == min(sds))
          if(length(indximproved) > 0)
          {
            Correctiontable[t+1,1+i] <- maxima[indximproved]
            maxdens[outlierindx] <- maxima[indximproved]
          }
        }
      }
      
      meanmaxdens <- mean(maxdens)
      ####plot density distribution before normalization
      
      if(meanmaxdens != 0) ##indicates that there were no datapoints available at all
      {
        cols <- c("black","red","blue","darkgreen","darkgrey","chocolate","burlywood4","cadetblue","chocolate4","chartreuse4","coral4","deeppink3") ### 12 colors for up to 12 replicates
        ind <- 0
        for(i in dataindices)
        {
          ind <- ind + 1
          temp <- data_list[[i]]
          temp <- subset(temp,qusm >= cfg_info$qusm_filter_norm & qupm >= cfg_info$qupm_filter_norm)
          temp <- temp[which(temp[,colindex[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
          if(norm_in_log2 == 1){density <- density(log2(temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]]))}else{density <- density(temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]])}
          
          if(ind == 1)
          {
            if(norm_in_log2 == 1)
            {
              plot(density,xlim=c(-2,2),ylim=c(0,densymax),col=cols[ind],main=paste("Rel.fcs.(log2) in ",cfg_info$condition_name_1," at ",cfg_info$temp_gradient[cfg_info$temp_order[t+1,g]]," C",sep=""))
            }else
            {
              plot(density,xlim=c(0,2),ylim=c(0,densymax),col=cols[ind],main=paste("Rel.fcs. in ",cfg_info$condition_name_1," at ",cfg_info$temp_gradient[cfg_info$temp_order[t+1,g]]," C",sep=""))
            }
          }else
          {
            lines(density,col=cols[ind])
          }
          ####indicate selected maxima
          if(norm_in_log2 == 1){maxima <- allmaxDensity(vals = log2(temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]]))}else{maxima <- allmaxDensity(vals = temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]])}
          for(m in maxima)
          {
            mindex <- which(density$x > m)[1]-1
            points(density$x[mindex],density$y[mindex],col=cols[ind])
          }
          segments(density$x[which(density$x>Correctiontable[t+1,1+i])[1]-1], density$y[which(density$x>Correctiontable[t+1,1+i])[1]-1], density$x[which(density$x>Correctiontable[t+1,1+i])[1]-1], par("usr")[3],lty=3, col = cols[ind])
        }
        legend("topright",col=cols[1:ind],legend = paste("Replicate ",1:length(dataindices),sep=""),lty=1)
        ####normalize Condition1 to mean of all replicates
        if(dontnormalize == F)
        {
          for(i in dataindices)
          {
            if(norm_in_log2 == 1){corrfactor <- 2^meanmaxdens / 2^Correctiontable[t+1,1+i]}else{corrfactor <- meanmaxdens / Correctiontable[t+1,1+i]}
            Correctiontable[t+1,1+cfg_info$num_samples+i] <- corrfactor
            normcolname <- column[which(dataindices == i)]
            
            data_list[[i]][,normcolname] <- data_list[[i]][,colindex[which(dataindices == i)]] * corrfactor
            
          }
        }else
        {
          for(i in dataindices)
          {
            corrfactor <- NA
            Correctiontable[t+1,1+cfg_info$num_samples+i] <- NA
            normcolname <- column[which(dataindices == i)]
            data_list[[i]][,normcolname] <- data_list[[i]][,colindex[which(dataindices == i)]]
          }
        }
        
        ####plot density distribution after normalization
        cols <- c("black","red","blue","darkgreen","darkgrey","chocolate","burlywood4","cadetblue","chocolate4","chartreuse4","coral4","deeppink3") ### 12 colors for up to 12 replicates
        ind <- 0
        for(i in dataindices)
        {
          temp <- data_list[[i]]
          temp <- subset(temp,qusm >= cfg_info$qusm_filter_norm & qupm >= cfg_info$qupm_filter_norm)
          temp <- temp[which(temp[,colindexnotnorm[which(dataindices == i)]] >= cfg_info$relfc_cutoff_norm),]
          ind <- ind + 1
          normcolname <- column[which(dataindices == i)]
          
          if(norm_in_log2 == 1){density <- density(log2(temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]]))}else{density <- density(temp[!is.na(temp[,colindex[which(dataindices == i)]]),colindex[which(dataindices == i)]])}
          
          if(ind == 1)
          {
            if(norm_in_log2 == 1)
            {
              plot(density,xlim=c(-2,2),ylim=c(0,densymax),col=cols[ind],main=paste("Normalized rel.fcs. (log2) in ", cfg_info$condition_name_1," at ",cfg_info$temp_gradient[cfg_info$temp_order[t+1,g]]," C",sep=""))
            }else
            {
              plot(density,xlim=c(0,2),ylim=c(0,densymax),col=cols[ind],main=paste("Normalized rel.fcs. in ", cfg_info$condition_name_1," at ",cfg_info$temp_gradient[cfg_info$temp_order[t+1,g]]," C",sep=""))
            }
          }else
          {
            lines(density,col=cols[ind])
          }
          segments(density$x[which(density$x>meanmaxdens)[1]-1], density$y[which(density$x>meanmaxdens)[1]-1], density$x[which(density$x>meanmaxdens)[1]-1], par("usr")[3],lty=3, col = cols[ind])
        }
        legend("topright",col=cols[1:ind],legend = paste("Replicate ",1:length(dataindices),sep=""),lty=1)
      }else #no available data points in the respective TMT channels thus skip normalization
      {
        for(i in dataindices)
        {
          corrfactor <- NA
          Correctiontable[t+1,1+cfg_info$num_samples+i] <- NA
          normcolname <- column[which(dataindices == i)]
          data_list[[i]][,normcolname] <- data_list[[i]][,colindex[which(dataindices == i)]]
        }
      }
      
    }
  }
  dev.off()
  
  output <- list()
  output[["data_list_norm"]] <- data_list
  output[["cond1_normalization_parameter_sheet"]] <- Correctiontable
  return(output)
}