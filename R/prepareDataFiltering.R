prepareDataFiltering = function(summary_data_combined,cfg_info,fordata = "both") ### prepare columns for data filtering
{
  ###Filter data:
  ###qupm 
  ###qusm
  ###relfc in condition1 or condition2 @ highest temp < 0.4
  ###R² of melting curve
  ###significance of protein based on < 50 % of temps in melting curve plateau <-- this filter will be applied on visualization (grey-out data points)
  ###significance of a protein was based on < 50% of temps at which the protein showed a high standard deviation <-- this filter will be applied on visualization (grey-out data points)
  ## determine in how many replicates qusmfilter is passed
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
    ## determine if protein was good quantified in at least two replicates per temperature gradient part: qusm
    npassedqusmfilter <- matrix(nrow = nrow(summary_data_combined[[ndata]]),ncol = cfg_info$num_gradientparts)
    for(i in 1:nrow(summary_data_combined[[ndata]]))
    {
      for(g in 1:cfg_info$num_gradientparts)
      {
        npassedqusmfilter[i,g] <- 0
        for(r in 1:cfg_info$replicates)
        {
          colindex <- which(grepl(paste("qusm.",r,".",g,sep=""),colnames(summary_data_combined[[ndata]])))
          if(!is.na(summary_data_combined[[ndata]][i,colindex]) && summary_data_combined[[ndata]][i,colindex] >= cfg_info$qusm_filter_plot){npassedqusmfilter[i,g] <- npassedqusmfilter[i,g] + 1}
        }
        if(npassedqusmfilter[i,g] >= 2){npassedqusmfilter[i,g] <- 1}else{npassedqusmfilter[i,g] <- 0} ###good quantification in at least two replicates
      }
    }
    npassedqusmfilter <- as.data.frame(npassedqusmfilter)
    colnames(npassedqusmfilter) <- c(paste("npassedqusmfilter_",1:cfg_info$num_gradientparts,sep=""))
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],npassedqusmfilter)
    
    ## determine if protein was good quantified in at least two replicates per temperature gradient part: qupm
    npassedqupmfilter <- matrix(nrow = nrow(summary_data_combined[[ndata]]),ncol = cfg_info$num_gradientparts)
    for(i in 1:nrow(summary_data_combined[[ndata]]))
    {
      for(g in 1:cfg_info$num_gradientparts)
      {
        npassedqupmfilter[i,g] <- 0
        for(r in 1:cfg_info$replicates)
        {
          colindex <- which(grepl(paste("qupm.",r,".",g,sep=""),colnames(summary_data_combined[[ndata]])))
          if(!is.na(summary_data_combined[[ndata]][i,colindex]) && summary_data_combined[[ndata]][i,colindex] >= cfg_info$qupm_filter_plot){npassedqupmfilter[i,g] <- npassedqupmfilter[i,g] + 1}
        }
        if(npassedqupmfilter[i,g] >= 2){npassedqupmfilter[i,g] <- 1}else{npassedqupmfilter[i,g] <- 0} ###good quantification in at least two replicates
      }
    }
    npassedqupmfilter <- as.data.frame(npassedqupmfilter)
    colnames(npassedqupmfilter) <- c(paste("npassedqupmfilter_",1:cfg_info$num_gradientparts,sep=""))
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],npassedqupmfilter)
    
    ## determine min mean relfc at highest temp for both conditions
    colcond1shighesttemp <- which(grepl(paste(cfg_info$temp_gradient[which(cfg_info$temp_gradient == as.character(max(as.numeric(cfg_info$temp_gradient))))],"C_cond1.",sep=""),colnames(summary_data_combined[[ndata]])))
    meanhighesttemp1 <- as.data.frame(rowMeans(summary_data_combined[[ndata]][,colcond1shighesttemp],na.rm=TRUE))
    colcond2shighesttemp <- which(grepl(paste(cfg_info$temp_gradient[which(cfg_info$temp_gradient == as.character(max(as.numeric(cfg_info$temp_gradient))))],"C_cond2.",sep=""),colnames(summary_data_combined[[ndata]])))
    meanhighesttemp2 <- as.data.frame(rowMeans(summary_data_combined[[ndata]][,colcond2shighesttemp],na.rm=TRUE))
    meanhighesttemp <- as.data.frame(matrix(ncol=1,nrow = nrow(meanhighesttemp1)))
    colnames(meanhighesttemp) <- "meanhightemp"
    for(i in 1:nrow(meanhighesttemp))
    {
      if(is.na(meanhighesttemp1[i,1]) && is.na(meanhighesttemp2[i,1]))
      {
        meanhighesttemp[i,1] <- NA
      }else
      {
        if(is.na(meanhighesttemp1[i,1]) && !is.na(meanhighesttemp2[i,1]))
        {
          meanhighesttemp[i,1] <- meanhighesttemp2[i,1]
        }else if(!is.na(meanhighesttemp1[i,1]) && is.na(meanhighesttemp2[i,1]))
        {
          meanhighesttemp[i,1] <- meanhighesttemp1[i,1]
        }else
        {
          if(meanhighesttemp1[i,1] <= meanhighesttemp2[i,1])
          {
            meanhighesttemp[i,1] <- meanhighesttemp1[i,1]
          }else
          {
            meanhighesttemp[i,1] <- meanhighesttemp2[i,1]
          }
        }
      }
    }
    meanhighesttemp[is.na(meanhighesttemp)] <- NA
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],meanhighesttemp)
    
    ## determine if significance is in majority based on temps where a plateau was reached
    SignificanceinPlateau <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
    colnames(SignificanceinPlateau) <- "Significance_in_plateau"
    for(i in 1:nrow(summary_data_combined[[ndata]]))
    {
      if(!is.na(summary_data_combined[[ndata]]$collation.ratio.start[i]) & !is.na(summary_data_combined[[ndata]]$collation.ratio.count[i]))
      {
        starttempindex <- which(as.numeric(cfg_info$temp_gradient) == cfg_info$temp_gradient[summary_data_combined[[ndata]]$collation.ratio.start[i]])
        endtempindex <- starttempindex+(summary_data_combined[[ndata]]$collation.ratio.count[i]-1)
        sigtemps <- as.numeric(cfg_info$temp_gradient[starttempindex:endtempindex])
        sigtempsinplateau <- which(sigtemps >= summary_data_combined[[ndata]]$Plateau[i])
        Plateaus <- summary_data_combined[[ndata]][i,which(grepl("Plateau_relfc_",colnames(summary_data_combined[[ndata]])))]
        Ratioplateaus <- log2(Plateaus[2]/Plateaus[1])
        if(!is.na(Ratioplateaus))
        {
          if((length(sigtempsinplateau)/length(sigtemps)) >= 0.5 & abs(Ratioplateaus) < log2(cfg_info$plateau_ratio_cutoff)) ###if more than 50 % of the sign. temps are sitting in the melting curve plateau AND the abundance difference between cond2 and cond1 plateaus < 2fold --> mark that protein as an outlier
          {
            SignificanceinPlateau[i,1] <- (length(sigtempsinplateau)/length(sigtemps))
          }
        }
      }
    }
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],SignificanceinPlateau)
    
    ## determine if proteins were found to be significantly affected but significance was based on >50% temps at which this protein has a too high standard deviation
    temp <- summary_data_combined[[ndata]][which(apply(summary_data_combined[[ndata]][,which(grepl("npassedqusmfilter",colnames(summary_data_combined[[ndata]]))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qusm filter
    temp <- temp[which(apply(temp[,which(grepl("npassedqupmfilter",colnames(temp))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qupm filter
      
    medianStDev <- colMedians(as.matrix(sapply(temp[,which(grepl("StDev_",colnames(temp)))], as.numeric)),na.rm=TRUE)
    StDevofStDev <- colSds(as.matrix(sapply(temp[,which(grepl("StDev_",colnames(temp)))], as.numeric)),na.rm=TRUE)
    factor <- ifelse(cfg_info$minstdev_factor_significance_function >= 3,cfg_info$minstdev_factor_significance_function,3)
    StDevcutoff <- medianStDev + factor*StDevofStDev
    ####now determine if significance was based on a temp at which this protein has a too high standard deviation
    Significance <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]]),""))
    colnames(Significance) <- "StDev_at_sig_temp_too_high"
    Significance$StDev_at_sig_temp_too_high <- as.character(Significance$StDev_at_sig_temp_too_high)
    for(i in 1:nrow(summary_data_combined[[ndata]]))
    {
      if(!is.na(summary_data_combined[[ndata]]$collation.ratio.start[i]) & !is.na(summary_data_combined[[ndata]]$collation.ratio.count[i]))
      {
        starttempindex <- which(as.numeric(cfg_info$temp_gradient) == cfg_info$temp_gradient[summary_data_combined[[ndata]]$collation.ratio.start[i]])
        endtempindex <- starttempindex+(summary_data_combined[[ndata]]$collation.ratio.count[i]-1)
        countoutlier <- 0
        outliertemps <- NULL
        for(sigtemp in starttempindex:endtempindex) ###determine if and on how many significant temps the stdev is higher than 3*median stdev of stdevs
        {
          if(!is.na(summary_data_combined[[ndata]][i,paste("StDev_",cfg_info$temp_gradient[sigtemp],"C",sep="")]))
          {
            if(summary_data_combined[[ndata]][i,paste("StDev_",cfg_info$temp_gradient[sigtemp],"C",sep="")] > StDevcutoff[sigtemp])
            {
              countoutlier <- countoutlier + 1
              outliertemps <- append(outliertemps,sigtemp)
            }
          }
        }
        
        ###if number of outlier stdevs > cutoff of num significant temps then remove protein from list
        if((countoutlier/summary_data_combined[[ndata]]$collation.ratio.count[i]) > cfg_info$relnum_sigtemps_stdev_outlier)
        {
          Significance$StDev_at_sig_temp_too_high[i] <- paste(cfg_info$temp_gradient[outliertemps],collapse = ",")
        }
      }
    }
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],Significance)
    
    ###calculate distance score for collated ratio and meanlog2ratioatreftemp
    ###for that filter data first
    temp <- summary_data_combined[[ndata]]
    temp <- temp[which(apply(temp[,which(grepl("npassedqusmfilter",colnames(temp))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qusm filter
    temp <- temp[which(apply(temp[,which(grepl("npassedqupmfilter",colnames(temp))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qupm filter
    temp <- subset(temp, temp[,"meanhightemp"] <= cfg_info$max_relfc_highesttemp_plot) ##relfc at highest temp
    temp <- subset(temp, rowMedians(as.matrix(temp[,which(grepl("Rsq_cond1",colnames(temp)))]),na.rm=T) >= cfg_info$Rsq_cutoff & rowMedians(as.matrix(temp[,which(grepl("Rsq_cond2",colnames(temp)))]),na.rm=T) >= cfg_info$Rsq_cutoff) ##R² filter
    temp <- temp[!duplicated(temp$Protein),]
    ###add protein of interest if they were removed by filters
    PoI <- as.matrix(subset(cfg_info$cfg_sheet$Proteins.of.interest,!is.na(cfg_info$cfg_sheet$Proteins.of.interest))) ### Proteins of interest 
    for(j in PoI)
    {
      if(length(which(temp$Protein == j)) == 0)
      {
        temp2 <- subset(summary_data_combined[[ndata]], Protein == j)
        temp <- rbind(temp, temp2)
      }
    }
    ###add column indicating for each protein if it fulfills all filtering criteria for stability change test
    summary_data_combined[[ndata]]$fullfil_all_criteria_collation <- ifelse(summary_data_combined[[ndata]]$Protein %in% temp$Protein,1,0)
    
    ###calculate distance score
    meancollatedratio <- mean(temp$collation.log2Ratio,na.rm=T)
    sdcollatedratio <- sd(temp$collation.log2Ratio,na.rm=T)
    summary_data_combined[[ndata]]$distancescore <- (summary_data_combined[[ndata]]$collation.log2Ratio-meancollatedratio)/sdcollatedratio
    
    ###BH-correction of pvalues for only those proteins which fullfil all filter
    temp$collation.pValue.corrected <- p.adjust(temp$collation.pValue,"BH",n = length(temp$collation.pValue))
    summary_data_combined[[ndata]] <- left_join(summary_data_combined[[ndata]],temp[,c("Protein","collation.pValue.corrected")],by=c("Protein"))
    
    ####get most ignificant pvalue at 37 °C for only those proteins which fullfil all filter criteria for 37 °C
    c <- which(grepl(paste("pValue_less",37,"C",sep=""),colnames(summary_data_combined[[ndata]])) | grepl(paste("pValue_greater",37,"C",sep=""),colnames(summary_data_combined[[ndata]])))
    columnname <- paste("pValue",37,"C",sep="")
    for(i in 1:nrow(summary_data_combined[[ndata]]))
    {
      summary_data_combined[[ndata]][i,columnname] <- min(summary_data_combined[[ndata]][i,c],na.rm=T)
      if(is.infinite(summary_data_combined[[ndata]][i,columnname])){summary_data_combined[[ndata]][i,columnname] <- NA}
    }
    temp <- summary_data_combined[[ndata]]
    temp <- temp[which(apply(temp[,which(grepl("npassedqusmfilter",colnames(temp))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qusm filter
    temp <- temp[which(apply(temp[,which(grepl("npassedqupmfilter",colnames(temp))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qupm filter
    temp <- temp[!duplicated(temp$Protein),]
    ###add protein of interest if they were removed by filters
    PoI <- as.matrix(subset(cfg_info$cfg_sheet$Proteins.of.interest,!is.na(cfg_info$cfg_sheet$Proteins.of.interest))) ### Proteins of interest 
    for(j in PoI)
    {
      if(length(which(temp$Protein == j)) == 0)
      {
        temp2 <- subset(summary_data_combined[[ndata]], Protein == j)
        temp <- rbind(temp, temp2)
      }
    }
    
    ###add column indicating for each protein if it fulfills all filtering criteria for abundance change test
    if(ndata == 1)
    {
      summary_data_combined[[ndata]]$fullfil_all_criteria_37 <- ifelse(summary_data_combined[[ndata]]$Protein %in% temp$Protein,1,0)
      
      meancollatedratio <- mean(temp$meanlog2Ratio37C,na.rm=T)
      sdcollatedratio <- sd(temp$meanlog2Ratio37C,na.rm=T)
      summary_data_combined[[ndata]]$distancescore37 <- (summary_data_combined[[ndata]]$meanlog2Ratio37C-meancollatedratio)/sdcollatedratio
      
      temp$pValue37C.corrected <- p.adjust(temp$pValue37C,"BH",n = length(temp$pValue37C))
      summary_data_combined[[ndata]] <- left_join(summary_data_combined[[ndata]],temp[,c("Protein","pValue37C.corrected")],by=c("Protein"))
    }

  }
  return(summary_data_combined)
}