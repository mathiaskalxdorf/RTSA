plotResults = function(summary_data_combined,norm_factors_meltingcurve=NA,cfg_info,fordata = "both",Internalizationcategories,Thermalshiftcategories,colors_melt_cond1,colors_melt_cond2)
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
  for(ndata in calcsetup)
  {
    default_par <- par() #save par
    
    Analysis <- summary_data_combined[[ndata]]
    
    for(t in cfg_info$temp_gradient)
    {
      c <- which(grepl(paste("pValue_less",t,"C",sep=""),colnames(Analysis)) | grepl(paste("pValue_greater",t,"C",sep=""),colnames(Analysis)))
      columnname <- paste("pValue",t,"C",sep="")
      for(i in 1:nrow(Analysis))
      {
        Analysis[i,columnname] <- min(Analysis[i,c],na.rm=T)
        if(is.infinite(Analysis[i,columnname])){Analysis[i,columnname] <- NA}
      }
    }
    
    if(ndata == 1)
    {
      pdf(file=paste("Result plots - with abundance effects.pdf",sep=""),useDingbats=FALSE)
    }
    if(ndata == 2)
    {
      pdf(file=paste("Result plots - without abundance effects.pdf",sep=""),useDingbats=FALSE)
    }
    
    ####save full analysis file
    Analysissave <- Analysis
    
    ####Plot mean StDev over all temps
    Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqusmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qusm filter
    Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqupmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qupm filter
    Analysis <- subset(Analysis, Analysis[,"meanhightemp"] <= cfg_info$max_relfc_highesttemp_plot) ##relfc at highest temp
    Analysis <- subset(Analysis, rowMedians(as.matrix(Analysis[,which(grepl("Rsq_cond1",colnames(Analysis)))]),na.rm=T) >= cfg_info$Rsq_cutoff & rowMedians(as.matrix(Analysis[,which(grepl("Rsq_cond2",colnames(Analysis)))]),na.rm=T) >= cfg_info$Rsq_cutoff) ##R² filter
    Analysis <- Analysis[!duplicated(Analysis$Protein),]
    ###add protein of interest if they were removed by filters
    PoI <- as.matrix(subset(cfg_info$cfg_sheet$Proteins.of.interest,!is.na(cfg_info$cfg_sheet$Proteins.of.interest))) ### Proteins of interest 
    for(j in PoI)
    {
      if(length(which(Analysis$Protein == j)) == 0)
      {
        temp2 <- subset(Analysissave, Protein == j)
        Analysis <- rbind(Analysis, temp2)
      }
    }
    
    meanSD <- as.matrix(Analysis[,which(grepl("meanStDev",colnames(Analysis)))])
    info <- boxplot(meanSD,main=paste("Standard deviation of log2(",cfg_info$condition_name_2,"/",cfg_info$condition_name_1,")",sep=""), xlab="StDev",ylim=c(0,2))
    text(x=1.3, labels = round(info$stats[3,1],digits=3), y=info$stats[3,1])
    mediansd <- info$stats[3,1]
    
    ###recalculation for distance score plotting
    
    ####Significance cut-off function all temps
    significancecutofffunction <- as.data.frame(matrix(ncol=2,nrow=100001))
    colnames(significancecutofffunction) <- c("x","y")
    meanlog2colratio <- mean(Analysis$collation.log2Ratio,na.rm=T)
    sdlog2colratio <- sd(Analysis$collation.log2Ratio,na.rm=T)
    significancecutofffunction$x <- seq(-30, 30, len = 100001)
    #borderminxleft <- ((-(cfg_info$minstdev_factor_significance_function*mediansd)-meanlog2colratio)/sdlog2colratio)
    #borderminxright <- (((cfg_info$minstdev_factor_significance_function*mediansd)-meanlog2colratio)/sdlog2colratio)
    #significancecutofffunction$y <- ifelse(significancecutofffunction$x < borderminxleft,(1/(abs(((sdlog2colratio*significancecutofffunction$x)+meanlog2colratio)+(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(significancecutofffunction$x > borderminxright,(1/(abs(((sdlog2colratio*significancecutofffunction$x)+meanlog2colratio)-(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),1000000))
    
    significancecutofffunction$y <- ifelse(significancecutofffunction$x < -(cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(significancecutofffunction$x+(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(significancecutofffunction$x > (cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(significancecutofffunction$x-(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),1000000))
    significancecutofffunction$x <- (significancecutofffunction$x-meanlog2colratio)/sdlog2colratio
    
    if(ndata == 1) ###prepare for abundance change plot
    {
      Analysis <- Analysissave
      Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqusmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qusm filter
      Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqupmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qupm filter
      
      ###add protein of interest if they were removed by filters
      PoI <- as.matrix(subset(cfg_info$cfg_sheet$Proteins.of.interest,!is.na(cfg_info$cfg_sheet$Proteins.of.interest))) ### Proteins of interest 
      for(j in PoI)
      {
        if(length(which(Analysis$Protein == j)) == 0)
        {
          temp2 <- subset(Analysissave, Protein == j)
          Analysis <- rbind(Analysis, temp2)
        }
      }
      
      
      #meanlog2colratio37 <- mean(Analysis$meanlog2Ratio37C,na.rm=T)
      #sdlog2colratio37 <- sd(Analysis$meanlog2Ratio37C,na.rm=T)
      
      meanSD <- as.matrix(Analysis[,which(grepl("StDev_37C",colnames(Analysis)))])
      info <- boxplot(meanSD,main=paste("37 \u00B0C standard deviation of log2(",cfg_info$condition_name_2,"/",cfg_info$condition_name_1,")",sep=""), xlab="StDev",ylim=c(0,2))
      text(x=1.3, labels = round(info$stats[3,1],digits=3), y=info$stats[3,1])
      mediansd37 <- info$stats[3,1]
      #mediansd37 <- (mediansd37-mean(summary_data_combined[[ndata]]$meanlog2Ratio37C,na.rm=T))/sd(summary_data_combined[[ndata]]$meanlog2Ratio37C,na.rm=T)
      
      ####Significance cut-off function 37 °C
      significancecutofffunction37 <- as.data.frame(matrix(ncol=2,nrow=100001))
      colnames(significancecutofffunction37) <- c("x","y")
      significancecutofffunction37$x <- seq(-30, 30, len = 100001)
      #borderminxleft37 <- ((-(cfg_info$minstdev_factor_significance_function*mediansd37)-meanlog2colratio37)/sdlog2colratio37)
      #borderminxright37 <- (((cfg_info$minstdev_factor_significance_function*mediansd37)-meanlog2colratio37)/sdlog2colratio37)
      #significancecutofffunction37$y <- ifelse(significancecutofffunction37$x < borderminxleft37,(0.1/(abs(((sdlog2colratio37*significancecutofffunction37$x)+meanlog2colratio37)+(cfg_info$minstdev_factor_significance_function*mediansd37))))+(-log10(cfg_info$pvalue_abundance_cutoff)),ifelse(significancecutofffunction37$x > borderminxright37,(0.1/(abs(((sdlog2colratio37*significancecutofffunction37$x)+meanlog2colratio37)-(cfg_info$minstdev_factor_significance_function*mediansd37))))+(-log10(cfg_info$pvalue_abundance_cutoff)),1000000))
      significancecutofffunction37$y <- ifelse(significancecutofffunction37$x < -(cfg_info$minstdev_factor_significance_function*mediansd37),(0.1/(abs(significancecutofffunction37$x+(cfg_info$minstdev_factor_significance_function*mediansd37))))+(-log10(cfg_info$pvalue_abundance_cutoff)),ifelse(significancecutofffunction37$x > (cfg_info$minstdev_factor_significance_function*mediansd37),(0.1/(abs(significancecutofffunction37$x-(cfg_info$minstdev_factor_significance_function*mediansd37))))+(-log10(cfg_info$pvalue_abundance_cutoff)),1000000))
      #significancecutofffunction37$x <- (significancecutofffunction37$x-meanlog2colratio37)/sdlog2colratio37
    }
    
    
    ##now filter data
    
    Analysis <- Analysissave
    Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqusmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qusm filter
    Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqupmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qupm filter
    Analysis <- subset(Analysis, Analysis[,"meanhightemp"] <= cfg_info$max_relfc_highesttemp_plot) ##relfc at highest temp
    Analysis <- subset(Analysis, rowMedians(as.matrix(Analysis[,which(grepl("Rsq_cond1",colnames(Analysis)))]),na.rm=T) >= cfg_info$Rsq_cutoff & rowMedians(as.matrix(Analysis[,which(grepl("Rsq_cond2",colnames(Analysis)))]),na.rm=T) >= cfg_info$Rsq_cutoff) ##R² filter ??? 3
    Analysis <- Analysis[!duplicated(Analysis$Protein),]
    
    ###add protein of interest if they were removed by filters
    PoI <- as.matrix(subset(cfg_info$cfg_sheet$Proteins.of.interest,!is.na(cfg_info$cfg_sheet$Proteins.of.interest))) ### Proteins of interest 
    for(j in PoI)
    {
      if(length(which(Analysis$Protein == j)) == 0)
      {
        temp <- subset(Analysissave, Protein == j)
        Analysis <- rbind(Analysis, temp)
      }
    }
    Analysis <- Analysis[order(Analysis$Protein),]
    
    ####Ranked MS1Intensities of all proteins
    MS1Intensities <- as.data.frame(matrix(nrow=nrow(Analysis),ncol=2))
    MS1Intensities[1]<- Analysis$Protein
    MS1Intensities[2]<- rowMeans(subset(Analysis, select = which(grepl("ms1intensity",colnames(Analysis)))), na.rm = TRUE)
    colnames(MS1Intensities) <- c("gene_name","ms1Intensities")
    MS1Intensities <- MS1Intensities[with(MS1Intensities, order(-ms1Intensities)), ]
    
    ####Calculate colors for boxes of significant proteins
    coloring <- matrix(nrow=nrow(Analysis),ncol=1) ###color for internalization
    coloring2 <- matrix(nrow=nrow(Analysis),ncol=1) ### color for thermal shift
    colnames(coloring) <- "Color_Internalization"
    colnames(coloring2) <- "Color_Thermalshift"
    for(inde in 1:nrow(Analysis)) #### go through all data and set color for both halfs of data point
    {
      #####determine values for internalization
      relfc37 <- 2^Analysis$meanlog2Ratio37C[inde]
      change37 <- round(relfc37 - 1,2) ###rel.abundance change to cond1
      
      if(!is.na(Analysis$pValue37C[inde]) && Analysis$pValue37C[inde] < cfg_info$pvalue_abundance_cutoff && abs(change37) >= cfg_info$abundance_cutoff)
      {
        if(abs(change37) < Internalizationcategories[1])
        {
          if(change37 > 0){coloring[inde] <- 1}
          if(change37 < 0){coloring[inde] <- -1}
        }else if(abs(change37) < Internalizationcategories[2])
        {
          if(change37 > 0){coloring[inde] <- 2}
          if(change37 < 0){coloring[inde] <- -2}
        }else if(abs(change37) >= Internalizationcategories[2])
        {
          if(change37 > 0){coloring[inde] <- 3}
          if(change37 < 0){coloring[inde] <- -3}
        }
      }else
      {
        coloring[inde] <- 0
      }
      
      #####determine values for thermal shift
      if(!is.na(Analysis$Mean_dTm[inde]) & !is.na(Analysis$dTm_pValue[inde]) & abs(Analysis$Mean_dTm[inde]) >= cfg_info$thermshift_cutoff & Analysis$dTm_pValue[inde] < cfg_info$pvalue_thermshift_cutoff) ###protein exceeding Thermalshiftcut and Thermalshiftpvaluecut
      {
        Meltshift <- Analysis$Mean_dTm[inde]
        if(abs(Meltshift) < Thermalshiftcategories[1])
        {
          if(Meltshift > 0){coloring2[inde] <- 1}
          if(Meltshift < 0){coloring2[inde] <- -1}
        }else if(abs(Meltshift) < Thermalshiftcategories[2])
        {
          if(Meltshift > 0){coloring2[inde] <- 2}
          if(Meltshift < 0){coloring2[inde] <- -2}
        }else if(abs(Meltshift) >= Thermalshiftcategories[2])
        {
          if(Meltshift > 0){coloring2[inde] <- 3}
          if(Meltshift < 0){coloring2[inde] <- -3}
        }
      }else
      {
        coloring2[inde] <- 0
      }
    }
    Analysis <- cbind(Analysis,coloring)
    Analysis <- cbind(Analysis,coloring2)
    
    if(ndata == 1) 
    {
      ####Draw volcano plot - all proteins
      par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
      with(Analysis, plot(distancescore, -log10(collation.pValue), col = "grey", pch=20,xlab = paste("Distance score (",cfg_info$condition_name_2,"/",cfg_info$condition_name_1,")",sep=""), ylab="Collated pValue,log10", main="Significantly affected proteins"))
      lines(ifelse(significancecutofffunction$x>=min(Analysis$distancescore,na.rm=T)-0.25 & significancecutofffunction$x<= max(Analysis$distancescore, na.rm=T)+0.25,significancecutofffunction$x,NA),ifelse(significancecutofffunction$y<=max(-log10(Analysis$collation.pValue),na.rm=T)+0.29,significancecutofffunction$y,NA),lty=2)
      par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)
      # Add vis box for significant prots
      significant <- subset(Analysis, -log10(collation.pValue) > ifelse(collation.log2Ratio < -(cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(collation.log2Ratio+(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(collation.log2Ratio > (cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(collation.log2Ratio-(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),1000))) 
      #significant <- subset(Analysis, -log10(collation.pValue) > ifelse(distancescore < borderminxleft,(1/abs(((sdlog2colratio*distancescore)+meanlog2colratio)+(cfg_info$minstdev_factor_significance_function*mediansd)))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(distancescore > borderminxright,(1/abs(((sdlog2colratio*distancescore)+meanlog2colratio)-(cfg_info$minstdev_factor_significance_function*mediansd)))+(-log10(cfg_info$min_pvalue_significance_function)),1000))) 
      
      if(cfg_info$plateau_stdev_plot_filter == 1)
      {
        significant2 <- subset(significant,StDev_at_sig_temp_too_high == "" & is.na(Significance_in_plateau))#Start_Collation_after_Meltpnt == 0 &  #remove proteins were collation started after melting point of protein or StDev at significant temp was too high
        ####color points based on which criteria data point was regarded as finally not significant
        if(nrow(significant) > 0)
        {
          for(inde in 1:nrow(significant)) 
          {
            # if(significant$Start_Collation_after_Meltpnt[inde] != 0)###if protein is not relevant as collation started after Mltpt
            # {
            #   with(significant[inde,], points(collation.log2Ratio,-log10(collation.pValue),col = "black",pch=19))
            # }
            if(significant$StDev_at_sig_temp_too_high[inde] != "")###if protein is not relevant as for >50% of sigtemps the StDev was too high
            {
              with(significant[inde,], points(distancescore,-log10(collation.pValue),col = "deeppink",pch=19,cex=0.7))
            }
            if(!is.na(significant$Significance_in_plateau[inde]))###if protein is not relevant as for >50% of sigtemps the StDev was too high
            {
              with(significant[inde,], points(distancescore,-log10(collation.pValue),col = "forestgreen",pch=7,cex=1))
            }
          }
        }
      }else
      {
        significant2 <- significant
      }
      cols1 <- colorRampPalette(c("#3A6C9A","white", "red"))(7) ###for internalization -3 to +3
      cols2 <- colorRampPalette(c("#3A6C9A","white", "red"))(7) ###for shift -3 to +3
      
      ###label significant
      if(nrow(significant2) > 0)
      {
        if(cfg_info$color_significant_by_shift_abundance == 1)
        {
          for(inde in 1:nrow(significant2)) ###only draw boxes for significant proteins fullfilling all critera
          {
            drawIntThermBox(significant2$distancescore[inde],-log10(significant2$collation.pValue[inde]),0.5,cols1[4+significant2$Color_Internalization[inde]],cols1[4+significant2$Color_Thermalshift[inde]],significant2$collation.pValue.corrected[inde] < cfg_info$min_pvalue_BH_corrected)
          }
        }else
        {
          points(significant2$distancescore,-log10(significant2$collation.pValue),col=ifelse(significant2$collation.pValue.corrected < cfg_info$min_pvalue_BH_corrected,"black","dimgray"),cex=1.5,pch=20)
          points(significant2$distancescore,-log10(significant2$collation.pValue),col=ifelse(is.na(significant2$collation.pValue.corrected),"dimgray",NA),cex=1.5,pch=20)
        }
      }
      
      # Label all points above cutoff line
      if(nrow(significant) > 0)
      {
        with(significant, textxy(distancescore, -log10(collation.pValue), labs=Protein, cex=0.5))
      }
      
      ###now add box and label for POI
      if(length(PoI) > 0 & cfg_info$color_significant_by_shift_abundance == 1)
      {
        for(prot in PoI)
        {
          inde <- which(Analysis$Protein == prot)
          if(length(inde) > 0)
          {
            if(!is.na(Analysis$distancescore[inde]) &!is.na(Analysis$collation.pValue[inde]))
            {
              if(cfg_info$color_significant_by_shift_abundance == 1)
              {
                if(prot %not in% significant2$Protein)
                {
                  drawIntThermBox(Analysis$distancescore[inde],-log10(Analysis$collation.pValue[inde]),0.5,cols1[4+Analysis$Color_Internalization[inde]],cols1[4+Analysis$Color_Thermalshift[inde]],Analysis$collation.pValue.corrected[inde] < cfg_info$min_pvalue_BH_corrected)
                }
              }else
              {
                points(significant2$distancescore,-log10(significant2$collation.pValue),col=ifelse(significant2$collation.pValue.corrected < cfg_info$min_pvalue_BH_corrected,"black","dimgray"),cex=1.5,pch=20)
              }
              textxy(Analysis$distancescore[inde], -log10(Analysis$collation.pValue[inde]), labs=Analysis$Protein[inde], cex=0.5,col="red")
            }
          }
          
        }
      }
      
      ###add legend
      if(cfg_info$color_significant_by_shift_abundance == 1)
      {
        drawLegendBox(cols1)
      }else
      {
        legend("topright",inset=c(-0.25,0),legend=c("BH-corrected","not corrected","not significant"),cex = 0.8,pt.cex=c(1.5,1.5,1) ,pch=c(20,20,20),col=c("black","dimgray","grey"),title = "Significance")
      }
      if(cfg_info$plateau_stdev_plot_filter == 1)
      {
        legend("bottomright",inset=c(-0.25,0),legend=c("Col. after MP","High StDev","Sig. in plateau"),cex = 0.8 ,pch=c(19,19,7),col=c("black","deeppink","forestgreen"),title = "Filter")
      }
      ####save list of significant proteins
      if(nrow(significant2) > 0)
      {
        saveExcel(Data = significant2,File = "Significant proteins - with abundance effect.xlsx")
      }
    }
    if(ndata == 2) 
    {
      #volcano plot with all proteins labeled
      par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
      with(Analysis, plot(distancescore, -log10(collation.pValue), col = "grey", pch=20,xlab = paste("Distance score (",cfg_info$condition_name_2,"/",cfg_info$condition_name_1,"),log2",sep=""), ylab="Collated pValue,log10", main="Significantly affected proteins"))
      lines(ifelse(significancecutofffunction$x>=min(Analysis$distancescore,na.rm=T)-0.25 & significancecutofffunction$x<= max(Analysis$distancescore, na.rm=T)+0.25,significancecutofffunction$x,NA),ifelse(significancecutofffunction$y<=max(-log10(Analysis$collation.pValue),na.rm=T)+0.29,significancecutofffunction$y,NA),lty=2)
      par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)
      # Add vis box for significant prots
      #significant <- subset(Analysis, -log10(collation.pValue) > ifelse(distancescore < borderminxleft,(1/abs(((sdlog2colratio*distancescore)+meanlog2colratio)+(cfg_info$minstdev_factor_significance_function*mediansd)))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(distancescore > borderminxright,(1/abs(((sdlog2colratio*distancescore)+meanlog2colratio)-(cfg_info$minstdev_factor_significance_function*mediansd)))+(-log10(cfg_info$min_pvalue_significance_function)),1000))) 
      significant <- subset(Analysis, -log10(collation.pValue) > ifelse(collation.log2Ratio < -(cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(collation.log2Ratio+(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(collation.log2Ratio > (cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(collation.log2Ratio-(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),1000))) 
      
      if(cfg_info$plateau_stdev_plot_filter == 1)
      {
        significant2 <- subset(significant,StDev_at_sig_temp_too_high == "" & is.na(Significance_in_plateau))#Start_Collation_after_Meltpnt == 0 &  #remove proteins were collation started after melting point of protein or StDev at significant temp was too high
        ####color points based on which criteria data point was regarded as finally not significant
        if(nrow(significant) > 0)
        {
          for(inde in 1:nrow(significant)) 
          {
            # if(significant$Start_Collation_after_Meltpnt[inde] != 0)###if protein is not relevant as collation started after Mltpt
            # {
            #   with(significant[inde,], points(collation.log2Ratio,-log10(collation.pValue),col = "black",pch=19))
            # }
            if(significant$StDev_at_sig_temp_too_high[inde] != "")###if protein is not relevant as for >50% of sigtemps the StDev was too high
            {
              with(significant[inde,], points(distancescore,-log10(collation.pValue),col = "deeppink",pch=19,cex=0.7))
            }
            if(!is.na(significant$Significance_in_plateau[inde]))###if protein is not relevant as for >50% of sigtemps the StDev was too high
            {
              with(significant[inde,], points(distancescore,-log10(collation.pValue),col = "forestgreen",pch=7,cex=1))
            }
          }
        }
      }else
      {
        significant2 <- significant
      }
      cols1 <- colorRampPalette(c("#3A6C9A","white", "red"))(7) ###for internalization -3 to +3
      cols2 <- colorRampPalette(c("#3A6C9A","white", "red"))(7) ###for shift -3 to +3
      if(nrow(significant2) > 0)
      {
        if(cfg_info$color_significant_by_shift_abundance == 1)
        {
          for(inde in 1:nrow(significant2)) ###only draw boxes for significant proteins fullfilling all critera
          {
            drawThermBox(significant2$distancescore[inde],-log10(significant2$collation.pValue[inde]),0.5,cols1[4+significant2$Color_Thermalshift[inde]],significant2$collation.pValue.corrected[inde] < cfg_info$min_pvalue_BH_corrected)
          }
        }else
        {
          points(significant2$distancescore,-log10(significant2$collation.pValue),col=ifelse(significant2$collation.pValue.corrected < cfg_info$min_pvalue_BH_corrected,"black","dimgray"),cex=1.5,pch=20)
          points(significant2$distancescore,-log10(significant2$collation.pValue),col=ifelse(is.na(significant2$collation.pValue.corrected),"dimgray",NA),cex=1.5,pch=20)
        }
      }
      
      # Label all points above cutoff line
      if(nrow(significant) > 0)
      {
        with(significant, textxy(distancescore, -log10(collation.pValue), labs=Protein, cex=0.5))
      }
      
      ###now add box and label for PoI
      if(length(PoI) > 0 & cfg_info$color_significant_by_shift_abundance == 1)
      {
        for(prot in PoI)
        {
          inde <- which(Analysis$Protein == prot)
          if(length(inde) > 0)
          {
            if(!is.na(Analysis$distancescore[inde]) &!is.na(Analysis$collation.pValue[inde]))
            {
              if(cfg_info$color_significant_by_shift_abundance == 1)
              {
                if(Analysis$Color_Thermalshift[inde] != 0 & prot %not in% significant2$Protein) ###add box only if thermal stability changed
                {
                  drawThermBox(Analysis$distancescore[inde],-log10(Analysis$collation.pValue[inde]),0.5,cols1[4+Analysis$Color_Thermalshift[inde]],Analysis$collation.pValue.corrected[inde] < cfg_info$min_pvalue_BH_corrected)
                }
              }else
              {
                if(Analysis$Color_Thermalshift[inde] != 0 & prot %not in% significant2$Protein) ###add box only if thermal stability changed
                {
                  drawThermBox(Analysis$distancescore[inde],-log10(Analysis$collation.pValue[inde]),0.5,cols1[4],Analysis$collation.pValue.corrected[inde] < cfg_info$min_pvalue_BH_corrected)
                }
              }
              textxy(Analysis$distancescore[inde], -log10(Analysis$collation.pValue[inde]), labs=Analysis$Protein[inde], cex=0.5,col="red")
            }
          }
          
        }
      }
      
      ###add legend
      if(cfg_info$color_significant_by_shift_abundance == 1)
      {
        drawLegendBox2(cols1)
      }else
      {
        legend("topright",inset=c(-0.25,0),legend=c("BH-corrected","not corrected","not significant"),cex = 0.8,pt.cex=c(1.5,1.5,1) ,pch=c(20,20,20),col=c("black","dimgray","grey"),title = "Significance")
      }
      if(cfg_info$plateau_stdev_plot_filter == 1)
      {
        legend("bottomright",inset=c(-0.25,0),legend=c("Sig. after MP","High StDev","Sig. in plateau","No therm. shift"),cex = 0.8 ,pch=c(19,19,7,20),col=c("black","deeppink","forestgreen","grey"),title = "Filter")
      }
      
      ####save list of significant proteins
      if(nrow(significant2) > 0)
      {
        saveExcel(Data = significant2,File = "Significant proteins - without abundance effect.xlsx")
      }
      
      
    }
    
    ####get significant proteins
    significant <- subset(Analysis, -log10(collation.pValue) > ifelse(collation.log2Ratio < -(cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(collation.log2Ratio+(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(collation.log2Ratio > (cfg_info$minstdev_factor_significance_function*mediansd),(1/(abs(collation.log2Ratio-(cfg_info$minstdev_factor_significance_function*mediansd))))+(-log10(cfg_info$min_pvalue_significance_function)),1000))) 
    #significant <- subset(Analysis, -log10(collation.pValue) > ifelse(distancescore < borderminxleft,(1/abs(((sdlog2colratio*distancescore)+meanlog2colratio)+(cfg_info$minstdev_factor_significance_function*mediansd)))+(-log10(cfg_info$min_pvalue_significance_function)),ifelse(distancescore > borderminxright,(1/abs(((sdlog2colratio*distancescore)+meanlog2colratio)-(cfg_info$minstdev_factor_significance_function*mediansd)))+(-log10(cfg_info$min_pvalue_significance_function)),1000))) 
    
    ###if proteins of interest were removed -> add them again
    for(j in PoI)
    {
      if(length(which(significant$Protein == j)) == 0)
      {
        temp <- subset(Analysis, Protein == j)
        significant <- rbind(significant, temp)
      }
    }
    significant <- significant[order(significant$Protein),]
    
    if(nrow(significant)>0)
    {
      if(cfg_info$normalize_melting_curve == 1)
      {
        ####Draw Mean-Melting curves for significant proteins with melting curve normalization
        drawMeltingcurves(cfg_info = cfg_info,data = significant,norm_factors_meltingcurve=norm_factors_meltingcurve,ndata=ndata,MS1Intensities = MS1Intensities,mediancurves = T,default_par = default_par)
        ####Draw melting curves for all replicates for significant proteins with melting curve normalization
        drawMeltingcurves(cfg_info = cfg_info,data = significant,norm_factors_meltingcurve=norm_factors_meltingcurve,ndata=ndata,MS1Intensities = MS1Intensities,mediancurves = F,default_par = default_par)
      }else
      {
        ####Draw Mean-Melting curves for significant proteins
        drawMeltingcurves(cfg_info = cfg_info,data = significant,ndata = ndata,MS1Intensities = MS1Intensities,mediancurves = T,default_par = default_par)
        ####Draw melting curves for all replicates for significant proteins
        drawMeltingcurves(cfg_info = cfg_info,data = significant,ndata = ndata,MS1Intensities = MS1Intensities,mediancurves = F,default_par = default_par)
        
      }
      
    }
    
    dev.off()
    
    if(ndata == 1) 
    {
      #use all good quantified proteins
      Analysis <- Analysissave
      Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqusmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qusm filter
      Analysis <- Analysis[which(apply(Analysis[,which(grepl("npassedqupmfilter",colnames(Analysis))),drop=F], MARGIN = 1, function(x) any(x == 1))),] #qupm filter
      
      ####Calculate colors for boxes of significant proteins
      coloring <- matrix(nrow=nrow(Analysis),ncol=1) ###color for internalization
      colnames(coloring) <- "Color_Internalization"
      for(inde in 1:nrow(Analysis)) #### go through all data and set color for both halfs of data point
      {
        #####determine values for internalization
        relfc37 <- 2^Analysis$meanlog2Ratio37C[inde]
        change37 <- relfc37 - 1 ###rel.abundance change to cond1
        
        if(!is.na(Analysis$pValue37C[inde]) && Analysis$pValue37C[inde] < cfg_info$pvalue_abundance_cutoff && abs(change37) >= cfg_info$abundance_cutoff)
        {
          if(abs(change37) < Internalizationcategories[1])
          {
            if(change37 > 0){coloring[inde] <- 1}
            if(change37 < 0){coloring[inde] <- -1}
          }else if(abs(change37) < Internalizationcategories[2])
          {
            if(change37 > 0){coloring[inde] <- 2}
            if(change37 < 0){coloring[inde] <- -2}
          }else if(abs(change37) >= Internalizationcategories[2])
          {
            if(change37 > 0){coloring[inde] <- 3}
            if(change37 < 0){coloring[inde] <- -3}
          }
        }else
        {
          coloring[inde] <- 0
        }
      }
      Analysis <- cbind(Analysis,coloring)
      
      ###Volcano plot showing abundance effects
      pdf("Induced abundance effect.pdf")
      par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=FALSE)
      with(Analysis, plot(meanlog2Ratio37C, -log10(pValue37C), col = "grey", pch=20,xlab = paste("37 \u00B0C mean ratio (",cfg_info$condition_name_2,"/",cfg_info$condition_name_1,"), log2",sep=""), ylab="37 \u00B0C pValue,log10", main="Significantly abundance regulated proteins"))
      lines(ifelse(significancecutofffunction37$x>=min(Analysis$meanlog2Ratio37C,na.rm=T)-0.25 & significancecutofffunction37$x<= max(Analysis$meanlog2Ratio37C, na.rm=T)+0.25,significancecutofffunction37$x,NA),ifelse(significancecutofffunction37$y<=max(-log10(Analysis$pValue37C),na.rm=T)+0.29,significancecutofffunction37$y,NA),lty=2)
      par(mar=c(5.1, 4.1, 4.1, 6.1), xpd=TRUE)
      # Add vis box for significant prots
      significant <- subset(Analysis, -log10(pValue37C) > ifelse(meanlog2Ratio37C < -(cfg_info$minstdev_factor_significance_function*mediansd37),(0.1/(abs(meanlog2Ratio37C+(cfg_info$minstdev_factor_significance_function*mediansd37))))+(-log10(cfg_info$pvalue_abundance_cutoff)),ifelse(meanlog2Ratio37C > (cfg_info$minstdev_factor_significance_function*mediansd37),(0.1/(abs(meanlog2Ratio37C-(cfg_info$minstdev_factor_significance_function*mediansd37))))+(-log10(cfg_info$pvalue_abundance_cutoff)),1000))) 
      cols1 <- colorRampPalette(c("#3A6C9A","white", "red"))(7) ###for internalization -3 to +3
      if(nrow(significant) > 0)
      {
        for(inde in 1:nrow(significant)) ###only draw boxes for significant proteins fullfilling all critera
        {
          drawThermBox(significant$meanlog2Ratio37C[inde],-log10(significant$pValue37C[inde]),0.5,cols1[4+significant$Color_Internalization[inde]],significant$pValue37C.corrected[inde] < cfg_info$min_pvalue_BH_corrected)
        }
        # Label all points above cutoff line
        with(significant, textxy(meanlog2Ratio37C, -log10(pValue37C), labs=Protein, cex=0.5))
        ####save list of significant proteins
        temp <- subset(significant,select=c("Protein","meanlog2Ratio37C","pValue37C"))
        temp <- cbind(temp,significant[,which(grepl("37C_",colnames(significant)))])
        temp <- cbind(temp,significant[,which(grepl("qusm\\.",colnames(significant)))])
        saveExcel(Data = temp,File = "Significant regulated proteins.xlsx")
      }
      ###add legend
      drawLegendBox3(cols1)
      
      dev.off()
    }
    
  }
}

drawIntThermBox <- function(x,y,cex,col1,col2,bold)##abundance and thermal shift
{
  plotminx <- par("usr")[1]
  plotmaxx <- par("usr")[2]
  plotminy <- par("usr")[3]
  plotmaxy <- par("usr")[4]
  sizex <- cex * (((plotmaxx-plotminx)*0.025))
  sizey <- cex * (((plotmaxy-plotminy)*0.025))
  if(!is.na(bold))
  {
    lwd <- ifelse(bold == T,2,1)
  }else
  {
    lwd <- 1
  }
  polygon(c(x-sizex, x, x, x-sizex,x-sizex), c(y-sizey, y-sizey, y+sizey, y+sizey, y-sizey), lty=1, col=col1,lwd=lwd)
  polygon(c(x,x+sizex,x+sizex,x,x), c(y-sizey, y-sizey, y+sizey, y+sizey, y-sizey), lty=1, col=col2,lwd=lwd)
}

drawThermBox <- function(x,y,cex,col1,bold)###thermal shift only
{
  plotminx <- par("usr")[1]
  plotmaxx <- par("usr")[2]
  plotminy <- par("usr")[3]
  plotmaxy <- par("usr")[4]
  if(!is.na(bold))
  {
    lwd <- ifelse(bold == T,2,1)
  }else
  {
    lwd <- 1
  }
  sizex <- cex * (((plotmaxx-plotminx)*0.025))
  sizey <- cex * (((plotmaxy-plotminy)*0.025))
  polygon(c(x-sizex, x+sizex, x+sizex, x-sizex,x-sizex), c(y-sizey, y-sizey, y+sizey, y+sizey, y-sizey), lty=1, col=col1,lwd=lwd)
}

drawLegendBox <- function(col)###abundance effect and thermal shift
{ 
  plotminx <- par("usr")[1]
  plotmaxx <- par("usr")[2]
  plotminy <- par("usr")[3]
  plotmaxy <- par("usr")[4]
  
  x <- plotmaxx
  y <- plotmaxy-(plotmaxy*0.15)
  
  startxbox <- plotmaxx + ((plotmaxx-plotminx)*0.08)
  sizexbox <- (((plotmaxx-plotminx)*0.03))
  
  startybox <- plotmaxy - 0.05*plotmaxy
  sizeybox <- (((plotmaxy-plotminy)*0.24))
  
  ###draw colors
  sizey <- sizeybox/7
  
  text(x+((((plotmaxx-plotminx)*0.24))*0.45),startybox + sizey,"Effect",cex=1.2)
  
  for(i in 1:7)
  {
    rect(startxbox,startybox-((i-1)*sizey),startxbox+(2*sizexbox),startybox-(i*sizey),col=col[8-i],border=NA)
    if(i == 1)
    {
      text(startxbox,startybox-((i-0.5)*sizey),"+50%",cex=0.8,pos=2,offset=0.1)
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),expression("">="2\u00B0C"),cex=0.8,pos=4,offset=0.1)
    }
    if(i == 2)
    {
      text(startxbox,startybox-((i-0.5)*sizey),"+25%",cex=0.8,pos=2,offset=0.1)
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"+1\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 3)
    {
      text(startxbox,startybox-((i-0.5)*sizey),"+10%",cex=0.8,pos=2,offset=0.1)
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"+0.5\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 5)
    {
      text(startxbox,startybox-((i-0.5)*sizey),"-10%",cex=0.8,pos=2,offset=0.1)
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"-0.5\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 6)
    {
      text(startxbox,startybox-((i-0.5)*sizey),"-25%",cex=0.8,pos=2,offset=0.1)
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"-1\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 7)
    {
      text(startxbox,startybox-((i-0.5)*sizey),"-50%",cex=0.8,pos=2,offset=0.1)
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),expression(""<="2\u00B0C"),cex=0.8,pos=4,offset=0.1)
    }
    
  }
  ###Draw box
  
  polygon(c(startxbox, startxbox+sizexbox, startxbox+sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizeybox, startybox-sizeybox, startybox), lty=1)
  polygon(c(startxbox+sizexbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox+sizexbox,startxbox+sizexbox), c(startybox, startybox, startybox-sizeybox, startybox-sizeybox, startybox), lty=1)
  
  text(startxbox+0.5*sizexbox,startybox-1.05*sizeybox,"Abundance",cex=0.8,pos=2,offset=0,srt = 90)
  text(startxbox+1.5*sizexbox,startybox-1.05*sizeybox,"Therm.shift",cex=0.8,pos=2,offset=0,srt = 90)
  
  ####draw legend for significance after BH correction
  
  startxbox <- plotmaxx + ((plotmaxx-plotminx)*0.02)
  startybox <- startybox-(15*sizey)
  
  text(x+((((plotmaxx-plotminx)*0.24))*0.45),startybox+sizey,"Significance",cex=1.2)
  
  polygon(c(startxbox, startxbox+sizexbox, startxbox+sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1,lwd=2)
  polygon(c(startxbox+sizexbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox+sizexbox,startxbox+sizexbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1,lwd=2)
  text(startxbox+(2*sizexbox),startybox-(0.5)*sizey,expression("with corr."),cex=0.8,pos=4,offset=0.1)
  
  startybox <- startybox-(1.5*sizey)
  polygon(c(startxbox, startxbox+sizexbox, startxbox+sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1)
  polygon(c(startxbox+sizexbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox+sizexbox,startxbox+sizexbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1)
  text(startxbox+(2*sizexbox),startybox-(0.5*sizey),expression("without corr."),cex=0.8,pos=4,offset=0.1)
  
  
  
}

drawLegendBox2 <- function(col)###only thermalshift
{ 
  plotminx <- par("usr")[1]
  plotmaxx <- par("usr")[2]
  plotminy <- par("usr")[3]
  plotmaxy <- par("usr")[4]
  
  x <- plotmaxx
  y <- plotmaxy-(plotmaxy*0.15)
  
  startxbox <- plotmaxx + ((plotmaxx-plotminx)*0.03)
  sizexbox <- (((plotmaxx-plotminx)*0.02))
  
  startybox <- plotmaxy - 0.05*plotmaxy
  sizeybox <- (((plotmaxy-plotminy)*0.24))
  
  ###draw colors
  sizey <- sizeybox/7
  
  for(i in 1:7)
  {
    rect(startxbox,startybox-((i-1)*sizey),startxbox+(2*sizexbox),startybox-(i*sizey),col=col[8-i],border=NA)
    if(i == 1)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),expression("">="2\u00B0C"),cex=0.8,pos=4,offset=0.1)
    }
    if(i == 2)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"+1\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 3)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"+0.5\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 5)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"-0.5\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 6)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"-1\u00B0C",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 7)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),expression(""<="2\u00B0C"),cex=0.8,pos=4,offset=0.1)
    }
    
  }
  ###Draw box
  polygon(c(startxbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizeybox, startybox-sizeybox, startybox), lty=1)
  
  text(startxbox+sizexbox,startybox-1.05*sizeybox,"Therm.shift",cex=0.8,pos=2,offset=0,srt = 90)
  
  ###draw legend for significance after BH correction
  
  startxbox <- plotmaxx + ((plotmaxx-plotminx)*0.02)
  startybox <- startybox-(15*sizey)
  
  text(x+((((plotmaxx-plotminx)*0.24))*0.45),startybox+sizey,"Significance",cex=1.2)
  
  polygon(c(startxbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1,lwd=2)
  text(startxbox+(2*sizexbox),startybox-(0.5)*sizey,expression("with corr."),cex=0.8,pos=4,offset=0.1)
  
  startybox <- startybox-(1.5*sizey)
  polygon(c(startxbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1)
  text(startxbox+(2*sizexbox),startybox-(0.5*sizey),expression("without corr."),cex=0.8,pos=4,offset=0.1)
  
}

drawLegendBox3 <- function(col)###only internalization
{ 
  plotminx <- par("usr")[1]
  plotmaxx <- par("usr")[2]
  plotminy <- par("usr")[3]
  plotmaxy <- par("usr")[4]
  
  x <- plotmaxx
  y <- plotmaxy-(plotmaxy*0.15)
  
  startxbox <- plotmaxx + ((plotmaxx-plotminx)*0.03)
  sizexbox <- (((plotmaxx-plotminx)*0.02))
  
  startybox <- plotmaxy - 0.05*plotmaxy
  sizeybox <- (((plotmaxy-plotminy)*0.24))
  
  ###draw colors
  sizey <- sizeybox/7
  
  text(x+((((plotmaxx-plotminx)*0.15))*0.45),startybox + sizey,"Effect",cex=1.2)
  
  for(i in 1:7)
  {
    rect(startxbox,startybox-((i-1)*sizey),startxbox+(2*sizexbox),startybox-(i*sizey),col=col[8-i],border=NA)
    if(i == 1)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),expression("">="+50%"),cex=0.8,pos=4,offset=0.1)
    }
    if(i == 2)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"+25%",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 3)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"+10%",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 5)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"-10%",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 6)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),"-25%",cex=0.8,pos=4,offset=0.1)
    }
    if(i == 7)
    {
      text(startxbox+(2*sizexbox),startybox-((i-0.5)*sizey),expression(""<="-50%"),cex=0.8,pos=4,offset=0.1)
    }
    
  }
  ###Draw box
  polygon(c(startxbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizeybox, startybox-sizeybox, startybox), lty=1)
  
  text(startxbox+sizexbox,startybox-1.05*sizeybox,"Abundance",cex=0.8,pos=2,offset=0,srt = 90)
  
  ###draw legend for significance after BH correction
  
  startxbox <- plotmaxx + ((plotmaxx-plotminx)*0.02)
  startybox <- startybox-(15*sizey)
  
  text(x+((((plotmaxx-plotminx)*0.24))*0.45),startybox+sizey,"Significance",cex=1.2)
  
  polygon(c(startxbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1,lwd=2)
  text(startxbox+(2*sizexbox),startybox-(0.5)*sizey,expression("with corr."),cex=0.8,pos=4,offset=0.1)
  
  startybox <- startybox-(1.5*sizey)
  polygon(c(startxbox, startxbox+2*sizexbox, startxbox+2*sizexbox, startxbox,startxbox), c(startybox, startybox, startybox-sizey, startybox-sizey, startybox), lty=1)
  text(startxbox+(2*sizexbox),startybox-(0.5*sizey),expression("without corr."),cex=0.8,pos=4,offset=0.1)
  
  
}

drawMeltingcurves <- function(cfg_info,data,norm_factors_meltingcurve=NA,ndata=NA,MS1Intensities,mediancurves=T,default_par)###data = data.frame containing relfc data, mediancurves defines if median of replicates is plotted or each replicate individually
{
  ylimmax <- matrix(nrow = nrow(data), ncol = 1)
  par(default_par)
  data <- data[order(data$Protein),]
  if(!is.na(norm_factors_meltingcurve) & !is.na(ndata)) ###normalize data based on melting curve normalization
  {
    normalized <- T
  }else
  {
    normalized <- F
  }
  for(i in 1:nrow(data))
  {
    par(fig=c(0,0.85,0,1), new=FALSE)
    startparam <- c("Pl"=0, "a"=550, "b"=10) 
    
    ####summarize relfcs for cur protein and cur replicate for all temps
    relfcscond1 <- matrix(ncol=length(cfg_info$temp_gradient),nrow=cfg_info$replicates*cfg_info$num_gradientparts)
    relfcscond2 <- relfcscond1
    for(t in cfg_info$temp_gradient)
    {
      cols <- which(grepl(paste(t,"C_cond1.",sep=""),colnames(data)))
      ind <- which(cfg_info$temp_gradient == t)
      if(length(cols) < cfg_info$replicates*cfg_info$num_gradientparts) ###only data for reftemp is available more than the number of replicates eg. 2 gradient parts 3 replicates -> 6 rows
      {
        temp <- c(as.numeric(data[i,cols]),rep(NA,(cfg_info$replicates*cfg_info$num_gradientparts)-length(cols)))
        relfcscond1[,ind] <- temp
      }else
      {
        relfcscond1[,ind] <- as.numeric(data[i,cols])
      }
      cols <- which(grepl(paste(t,"C_cond2.",sep=""),colnames(data)))
      if(length(cols) < cfg_info$replicates*cfg_info$num_gradientparts) ###only data for reftemp is available more than the number of replicates eg. 2 gradient parts 3 replicates -> 6 rows
      {
        temp <- c(as.numeric(data[i,cols]),rep(NA,(cfg_info$replicates*cfg_info$num_gradientparts)-length(cols)))
        relfcscond2[,ind] <- temp
      }else
      {
        relfcscond2[,ind] <- as.numeric(data[i,cols])
      }
    }
    
    if(normalized == T) ###normalize data based on melting curve normalization
    {
      if(ndata == 1)
      {
        relfcscond1[1:cfg_info$replicates,] <- relfcscond1[1:cfg_info$replicates,]*as.matrix(norm_factors_meltingcurve[[1]])
        relfcscond2[1:cfg_info$replicates,] <- relfcscond2[1:cfg_info$replicates,]*as.matrix(norm_factors_meltingcurve[[2]])
      }
      if(ndata == 2)
      {
        relfcscond1[1:cfg_info$replicates,] <- relfcscond1[1:cfg_info$replicates,]*as.matrix(norm_factors_meltingcurve[[3]])
        relfcscond2[1:cfg_info$replicates,] <- relfcscond2[1:cfg_info$replicates,]*as.matrix(norm_factors_meltingcurve[[4]])
      }
    }
    
    if(mediancurves == T)
    {
      ###calculate StDev for individual temp for errorbars in plots
      StDevcond1 <- matrix(,nrow = 1,ncol = length(cfg_info$temp_gradient))
      StDevcond2 <- matrix(,nrow = 1,ncol = length(cfg_info$temp_gradient))
      ###Condition1 and Condition2
      for(t in cfg_info$temp_gradient)
      {
        ind <- which(cfg_info$temp_gradient == t)
        StDevcond1[ind] <- 1.96*sd(relfcscond1[,ind],na.rm = TRUE)/sqrt(cfg_info$replicates) #1.96 -> 95% confedence interval
        StDevcond2[ind] <- 1.96*sd(relfcscond2[,ind],na.rm = TRUE)/sqrt(cfg_info$replicates) #1.96 -> 95% confedence interval
      }
      ####calculate median at each temp
      relfcscond1 <- as.data.frame(t(colMedians(relfcscond1,na.rm=T)))
      relfcscond2 <- as.data.frame(t(colMedians(relfcscond2,na.rm=T)))
      
    }
    ####determine ymax for cond2
    
    ylimmax[i,1] <- max(append(as.numeric(relfcscond1),as.numeric(relfcscond2)),na.rm=T)
    if(is.na(ylimmax[i]) || ylimmax[i] <= 1.2){ylimmax[i] <- 1.2}else{ylimmax[i] <- ylimmax[i] + 0.2}
    
    ####plot melting curves
    x <- as.numeric(cfg_info$temp_gradient)
    
    ###Preparation for Meltingpoints
    if(mediancurves == T)
    {
      Meltpoints <- as.data.frame(matrix(ncol=2,nrow=1))
      Rsquare <- as.data.frame(matrix(ncol=2,nrow=1))
      r <- 1
    }else
    {
      Meltpoints <- as.data.frame(matrix(ncol=2,nrow=cfg_info$replicates))
      Rsquare <- as.data.frame(matrix(ncol=2,nrow=cfg_info$replicates))
      r <- cfg_info$replicates
    }
    colnames(Meltpoints) <- c("cond1","cond2")
    colnames(Rsquare) <- c("cond1","cond2")
    
    for(j in 1:r)
    {
      ###plot meltcurve of condition1
      y <- as.numeric(relfcscond1[j,])
      if(j == 1)
      {
        if(normalized == T){main <- paste("Normalized melting curves of ",data$Protein[i],sep="")}else{main <- paste("Melting curves of ",data$Protein[i],sep="")}
        
        plot(x,y, col = colors_melt_cond1[j],pch=20,main=main,ylim=c(0,ylimmax[i]),xlab = "Temperature [\u00B0C]",ylab = "Fraction non-denatured")
        mtext(paste("col.ratio: ",round(data$collation.log2Ratio[i],digits=1)," col.pval: ",round(-log10(data$collation.pValue[i]),digits=1),sep=""))
        if(length(relfcscond1[,1]) > r) #### reference temp can contain more replicates - plot these extra data points
        {
          for(extra in 1:(cfg_info$num_gradientparts-1))
          {
            points(x[1],relfcscond1[r*extra+j,1], col = colors_melt_cond1[j],pch=20)
          }
        }
      }else
      {
        points(x,y, col = colors_melt_cond1[j],pch=20)
        if(length(relfcscond1[,1]) > r) #### reference temp can contain more replicates - plot these extra data points
        {
          for(extra in 1:(cfg_info$num_gradientparts-1))
          {
            points(x[1],relfcscond1[r*extra+j,1], col = colors_melt_cond1[j],pch=20)
          }
        }
      }
      if(mediancurves == T){errbar(x,y, y+StDevcond1, y-StDevcond1, add=T, pch=1, cap=.01)}
      plotminx <- par("usr")[1]
      plotmaxx <- par("usr")[2]
      plotminy <- par("usr")[3]
      plotmaxy <- par("usr")[4]
      y[1] <- median(relfcscond1[,1],na.rm=T)
      fit <- fitSigmoidTR2(x,y,startparam,10,y[1])
      new = data.frame(x = seq(min(x),max(x),len=500))
      if(class(fit) != "try-error")
      {
        lines(new$x,predict(fit,newdata=new),col = colors_melt_cond1[j],lwd=2)
        inflec <- 0.5*y[1]
        if(ndata == 1)###data with abundance effect require a different melting point determination
        {
          inflec <- y[1] - ((y[1] - predict(fit,data.frame(x=100)))/2)
          try(Meltpoints$cond1[j] <- (new$x[which(predict(fit,newdata=new)<=inflec)[1]]+new$x[which(predict(fit,newdata=new)<=inflec)[1]-1])/2,silent=TRUE)
        }else
        {
          pars <- coefficients(fit)
          a <- pars[["a"]]
          b <- pars[["b"]]
          pl <- pars[["Pl"]]
          Meltpoints$cond1[j] <- a / (b - log((1-pl)/(1/2 - pl) - 1))
        }
        
        if(!is.na(Meltpoints$cond1[j]))
        {
          #### indicate melting point ####
          segments(plotminx, inflec, Meltpoints$cond1[j], inflec,lty=2, col = colors_melt_cond1[j])
          segments(Meltpoints$cond1[j], inflec, Meltpoints$cond1[j], plotminy,lty=2, col = colors_melt_cond1[j])
        }
        ###calculate R²
        distres <- matrix(nrow = 1, ncol = length(cfg_info$temp_gradient))
        for(p in 1:(length(cfg_info$temp_gradient)-1))
        {
          n <- 1
          distres[n,p] <- ((y[p]-((predict(fit,newdata=new)[which(new$x>as.numeric(cfg_info$temp_gradient[p]))[1]]+predict(fit,newdata=new)[which(new$x>as.numeric(cfg_info$temp_gradient[p]))[1]-1])/2)))   
        }
        p <- length(cfg_info$temp_gradient)
        n <- 1
        distres[n,p] <- ((y[p]-(predict(fit,newdata=new)[500])))   
        disttot <- matrix(nrow = 1, ncol = length(cfg_info$temp_gradient))
        for(p in 1:length(cfg_info$temp_gradient))
        {
          n <- 1
          disttot[n,p] <- (y[p]-mean(y,na.rm=TRUE))
        }
        Rsquare$cond1[j] <- (1-(sum(distres[1,]^2,na.rm=T)/sum(disttot[1,]^2,na.rm=T)))
      }
      ###plot meltcurve of cond2
      y <- as.numeric(relfcscond2[j,])
      
      points(x,y, col = colors_melt_cond2[j],pch=20)
      if(length(relfcscond2[,1]) > r) #### reference temp can contain more replicates - plot these extra data points
      {
        for(extra in 1:(cfg_info$num_gradientparts-1))
        {
          points(x[1],relfcscond2[r*extra+j,1], col = colors_melt_cond2[j],pch=20)
        }
      }
      if(mediancurves == T){errbar(x,y, y+StDevcond2, y-StDevcond2, add=T, pch=1, cap=.01)}
      # ###check for significant abundance change at reftemp
      # if(data$pValue37C[i] > cfg_info$pvalue_abundance_cutoff | abs(round(2^data$meanlog2Ratio37C[i]-1,2)) < cfg_info$abundance_cutoff)##no significant abundance change
      # {
      #   y[1] <- 1
      # }else
      # {
      #   y[1] <- median(relfcscond2[,1],na.rm=T)
      # }
      y[1] <- median(relfcscond2[,1],na.rm=T)
      fit <- fitSigmoidTR2(x,y,startparam,10,y[1])
      new = data.frame(x = seq(min(x),max(x),len=500))
      if(class(fit) != "try-error")
      {
        lines(new$x,predict(fit,newdata=new),col = colors_melt_cond2[j],lwd=2)
        inflec <- 0.5*y[1]
        if(ndata == 1)###data with abundance effect require a different melting point determination
        {
          inflec <- y[1] - ((y[1] - predict(fit,data.frame(x=100)))/2)
          try(Meltpoints$cond2[j] <- (new$x[which(predict(fit,newdata=new)<=inflec)[1]]+new$x[which(predict(fit,newdata=new)<=inflec)[1]-1])/2,silent=TRUE)
        }else
        {
          pars <- coefficients(fit)
          a <- pars[["a"]]
          b <- pars[["b"]]
          pl <- pars[["Pl"]]
          Meltpoints$cond2[j] <- a / (b - log((1-pl)/(1/2 - pl) - 1))
        }
        
        if(!is.na(Meltpoints$cond2[j]))
        {
          #### indicate melting point ####
          segments(plotminx, inflec, Meltpoints$cond2[j], inflec,lty=2, col = colors_melt_cond2[j])
          segments(Meltpoints$cond2[j], inflec, Meltpoints$cond2[j], plotminy,lty=2, col = colors_melt_cond2[j])
        }
        ###calculate R²
        distres <- matrix(nrow = 1, ncol = length(cfg_info$temp_gradient))
        for(p in 1:(length(cfg_info$temp_gradient)-1))
        {
          n <- 1
          distres[n,p] <- ((y[p]-((predict(fit,newdata=new)[which(new$x>as.numeric(cfg_info$temp_gradient[p]))[1]]+predict(fit,newdata=new)[which(new$x>as.numeric(cfg_info$temp_gradient[p]))[1]-1])/2)))   
        }
        p <- length(cfg_info$temp_gradient)
        n <- 1
        distres[n,p] <- ((y[p]-(predict(fit,newdata=new)[500])))   
        disttot <- matrix(nrow = 1, ncol = length(cfg_info$temp_gradient))
        for(p in 1:length(cfg_info$temp_gradient))
        {
          n <- 1
          disttot[n,p] <- (y[p]-mean(y,na.rm=TRUE))
        }
        Rsquare$cond2[j] <- (1-(sum(distres[1,]^2,na.rm=T)/sum(disttot[1,]^2,na.rm=T)))
      }
    }
    
    #legend
    legend("topright", lty=1,lwd=2,col=c(colors_melt_cond1[1:r],colors_melt_cond2[1:r]), c(paste(cfg_info$condition_name_1," ",round(Meltpoints$cond1,digits=1),"\u00B0C (R\u00B2 ",round(Rsquare$cond1,digits=2),")",sep=""), paste(cfg_info$condition_name_2," ",round(Meltpoints$cond2,digits=1),"\u00B0C (R\u00B2 ",round(Rsquare$cond2,digits=2),")",sep="")), bty="o",cex=0.8, box.col="black",title = "Samples and Melt.pts.")
    
    ##### visualize thermal-shift#####
    if(mediancurves == T)
    {
      if(!is.na(Meltpoints$cond1) && !is.na(Meltpoints$cond2))
      {
        segments(Meltpoints$cond1, 0.03, Meltpoints$cond2, 0.03,lty=2, col = "dimgrey")
        if(!is.na(data$Mean_dTm[i]) & !is.na(data$dTm_pValue[i]) & abs(data$Mean_dTm[i]) >= cfg_info$thermshift_cutoff & data$dTm_pValue[i] < cfg_info$pvalue_thermshift_cutoff) ###protein exceeding Thermalshiftcut and Thermalshiftpvaluecut
        {
          textxy(Meltpoints$cond1 + (data$Mean_dTm[i]/2), 0.06, labs=paste("dif ",round(data$Mean_dTm[i],digits=1),sep=""), cex=0.6,offset = 0, col = "black" )
          textxy(Meltpoints$cond1 + (data$Mean_dTm[i]/2), 0.0, labs=paste("-log(p) ",round(-log10(data$dTm_pValue[i]),digits=1),sep=""), cex=0.6,offset = 0, col = "black" )
        }
      }
      ##### visualize internalization#####
      if(!is.na(relfcscond2[1,1]))
      {
        segments(x[1], 1, x[1], relfcscond2[1,1],lty=2, col = "dimgrey")
        
        if(!is.na(data[i,paste("pValue",cfg_info$temp_gradient[1],"C",sep="")]) && -log10(data[i,paste("pValue",cfg_info$temp_gradient[1],"C",sep="")]) >= -log10(cfg_info$pvalue_abundance_cutoff) && abs(relfcscond2[1,1]-1) >= cfg_info$abundance_cutoff) ##draw only if significant internalization
        {
          textxy(x[1]+2.6, 1 - ((1-relfcscond2[1,1])/2), labs=paste("relfc ",round((relfcscond2[1,1]-1)*100,digits=0),"%\n -log(p) ",round(-log10(data[i,paste("pValue",cfg_info$temp_gradient[1],"C",sep="")]),digits=1),sep=""), cex=0.6,offset = 0, col = "black" )   
        }
      }
    }
    
    #### highlight which datapoints were collated depending on selected collation function
    
    if(cfg_info$pVal_collation_function %in% c(0,2,3)) ###collate all available data
    {
      start <- data[i,"collation.ratio.start"]
      if(length(start) > 0)
      {
        if(!is.na(start))
        {
          for(count in 1:data[i,"collation.ratio.count"])
          {
            index <- start + count - 1
            rect(x[index]-1,min(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)-0.15,x[index]+1,max(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)+0.15,lty=3)
          }
        }
      }
    }
    
    if(cfg_info$pVal_collation_function == 1) ###collation with sliding-window approach
    {
      #### highlight which datapoints were collated for ratios
      start <- data[i,"collation.ratio.start"]
      if(!is.na(start) > 0 & !is.na(data[i,"collation.ratio.count"]))
      {
        if(!is.na(start))
        {
          for(count in 1:data[i,"collation.ratio.count"])
          {
            index <- start + count - 1
            rect(x[index]-1,min(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)-0.15,x[index]+1,max(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)+0.15,lty=3)
          }
        }
      }
    
      #### highlight which datapoints were collated for pvalues
      start <- data[i,"collation.pValue.start"]
      if(length(start)>0)
      {
        if(!is.na(start))
        {
          for(count in 1:data[i,"collation.pvalue.count"])
          {
            index <- start + count - 1
            rect(x[index]-1,min(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)-0.15,x[index]+1,max(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)+0.15,lty=3,border="red")
          }
        }
      }
    }
    
    # if(cfg_info$pVal_collation_function == 0 | cfg_info$pVal_collation_function == 3) ###collate all significant (Browns method or PECA)
    # {
    #   start <- which(cfg_info$temp_gradient == data[i,"collation.ratio.start"])
    #   if(length(start) > 0)
    #   {
    #     for(count in 1:data[i,"collation.ratio.count"])
    #     {
    #       index <- start + count - 1
    #       rect(x[index]-1,min(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)-0.15,x[index]+1,max(append(as.numeric(relfcscond1[,index]),as.numeric(relfcscond2[,index])),na.rm=T)+0.15,lty=3)
    #     }
    #   }
    # }
    
    ####add MS1Intensity plot
    par(fig=c(0.7,1,0.4,1),new=TRUE)
    plot(MS1Intensities$ms1Intensities,cex.main=0.8,  col="royalblue4",type = 'p', pch = 16,xlab='', ylab='', cex.axis=0.5,yaxt='n',xaxt='n')
    title("ms1Intensity", line = 0.7,cex.main = 0.8)
    axis(4, at = 1:10, label = rep("", 10), tck = -0.05)
    axis(4, at = 1:10, line = -0.9, lwd = 0, cex.axis = 0.9)
    mtext("ms1Intensity", 4, line = 0.9) 
    if(!is.na(MS1Intensities$ms1Intensities[which(MS1Intensities$gene_name==data$Protein[i])]))
    {
      points(which(MS1Intensities$gene_name==data$Protein[i]),MS1Intensities$ms1Intensities[which(MS1Intensities$gene_name==data$Protein[i])],col="red",pch=20,cex=2)
      par(xpd=T)
      textxy(which(MS1Intensities$gene_name==data$Protein[i]),MS1Intensities$ms1Intensities[which(MS1Intensities$gene_name==data$Protein[i])],data$Protein[i])
    }
    par(xpd=F)
    ####add qusm plot
    par(fig=c(0.7,1,0.0,0.6),new=TRUE)
    barplotdata <- as.data.frame(matrix(nrow = cfg_info$replicates*cfg_info$num_gradientparts,ncol=1))
    for(g in 1:cfg_info$num_gradientparts)
    {
      for(r in 1:cfg_info$replicates)
      {
        barplotdata[((g-1)*cfg_info$replicates)+r,1] <- data[i,which(grepl(paste("qusm.",r,".",g,sep=""),colnames(data)))]
        rownames(barplotdata)[((g-1)*cfg_info$replicates)+r] <- paste("Rep",r,".",g,sep="")
      }
    }
    colnames(barplotdata) <- "data"
    mids <- barplot(barplotdata$data,ylab='', cex.axis=0.1,yaxt='n',las=2)
    axis(1, at=mids, labels=rownames(barplotdata), las=3,cex.axis=0.5)
    axis(4, at = 0:100, label = rep("", 101), tck = -0.05)
    axis(4, at = 0:100, line = -0.9, lwd = 0, cex.axis = 0.9)
    mtext("qusm", 4, line = 0.9)
  }
}