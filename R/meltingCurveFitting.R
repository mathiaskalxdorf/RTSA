meltingCurveFitting = function(summary_data_combined,cfg_info,fordata = "both",progressbar = T) ###fit melting curves, calculate melting points and rsquares as well as determine melting plateaus
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
  norm_factors_meltingcurve_output <- list()
  
  for(ndata in calcsetup)
  {
    temps <- as.numeric(cfg_info$temp_gradient)
    startparam <- c("Pl"=0, "a"=550, "b"=10)
    Meltingpoints <- as.data.frame(matrix(ncol=2*length(unique(cfg_info$dataset_to_gradientpart$Replicate)),nrow=nrow(summary_data_combined[[ndata]])))
    RSqares <- Meltingpoints
    colnames(Meltingpoints) <- c(paste("Meltpnt_cond1_",unique(cfg_info$dataset_to_gradientpart$Replicate),sep=""),paste("Meltpnt_cond2_",unique(cfg_info$dataset_to_gradientpart$Replicate),sep=""))
    colnames(RSqares) <- c(paste("Rsq_cond1_",unique(cfg_info$dataset_to_gradientpart$Replicate),sep=""),paste("Rsq_cond2_",unique(cfg_info$dataset_to_gradientpart$Replicate),sep=""))
    numreplicates <- length(unique(cfg_info$dataset_to_gradientpart$Replicate))
    max <- nrow(summary_data_combined[[ndata]])*numreplicates
    counter <- 0
    
    if(cfg_info$normalize_melting_curve == 1)
    {
      ####additional data normalization for melting curve fitting comparable to original TPP publication
      ###determine normalizationfactors per temperature to achieve optimal melting curve fits
      
      ###select those proteins which were good quantified in all gradientparts and go below 0.2 at highest temp
      qusmcols <- which(grepl("^qusm",colnames(summary_data_combined[[ndata]])))
      select <- NULL
      colcond1shighesttemp <- which(grepl(paste(cfg_info$temp_gradient[which(cfg_info$temp_gradient == as.character(max(as.numeric(cfg_info$temp_gradient))))],"C_cond1.",sep=""),colnames(summary_data_combined[[ndata]])))
      colcond2shighesttemp <- which(grepl(paste(cfg_info$temp_gradient[which(cfg_info$temp_gradient == as.character(max(as.numeric(cfg_info$temp_gradient))))],"C_cond2.",sep=""),colnames(summary_data_combined[[ndata]])))
      for(i in 1:nrow(summary_data_combined[[ndata]]))
      {
        if(any(is.na(summary_data_combined[[ndata]][i,qusmcols])) == F)
        {
          if(any(summary_data_combined[[ndata]][i,qusmcols] < cfg_info$qusm_filter_norm) == F)
          {
            meanhighesttemp1 <- mean(as.numeric(summary_data_combined[[ndata]][i,colcond1shighesttemp]),na.rm=TRUE)
            meanhighesttemp2 <- mean(as.numeric(summary_data_combined[[ndata]][i,colcond2shighesttemp]),na.rm=TRUE)
            if(!is.na(meanhighesttemp1) & !is.na(meanhighesttemp2))
            {
              if(meanhighesttemp1 < 0.2 & meanhighesttemp2 < 0.2)
              {
                select <- append(select,i)
              }
            }
          }
        }
      }
      normProts <- summary_data_combined[[ndata]][select,]
      
      ####now generate melting curve fit based on all data
      
      relfcsctrl <- matrix(ncol=length(cfg_info$temp_gradient),nrow=0)
      relfcstreat <- matrix(ncol=length(cfg_info$temp_gradient),nrow=0)
      temprelfcctrl <- matrix(ncol=length(cfg_info$temp_gradient),nrow=1)
      temprelfctreat <- matrix(ncol=length(cfg_info$temp_gradient),nrow=1)
      for(i in 1:nrow(normProts))
      {
        for(r in unique(cfg_info$dataset_to_gradientpart$Replicate)) ###for each replicate
        {
          for(t in cfg_info$temp_gradient)
          {
            cols <- which(grepl(paste(t,"C_cond1.",r,sep=""),colnames(normProts)))
            ind <- which(cfg_info$temp_gradient == t)
            temprelfcctrl[1,ind] <- median(as.numeric(normProts[i,cols]),na.rm=T)
            cols <- which(grepl(paste(t,"C_cond2.",r,sep=""),colnames(normProts)))
            temprelfctreat[1,ind] <- median(as.numeric(normProts[i,cols]),na.rm=T)
          }
          relfcsctrl <- rbind(relfcsctrl,temprelfcctrl)
          relfcstreat <- rbind(relfcstreat,temprelfctreat)
        }
      }
      relfcstreat[,1] <- 1
      relfcs <- rbind(relfcsctrl,relfcstreat)
      ###median melting curve data from cond1 and cond2
      relfcsmedian <- colMedians(relfcs,na.rm=T)
      relfcsmedianctrl <- colMedians(relfcsctrl,na.rm=T)
      relfcsmediantreat <- colMedians(relfcstreat,na.rm=T)
      ###fit median melting curve
      fit <- fitSigmoidTR(temps,relfcsmedian,startparam,10,fixT0=TRUE)
      new = data.frame(x = seq(min(temps),max(temps),len=500))
      
      norm_factors_meltingcurve_ctrl <- matrix(nrow=length(unique(cfg_info$dataset_to_gradientpart$Replicate)),ncol=length(cfg_info$temp_gradient)) #normalizationfactors per replicate
      rownames(norm_factors_meltingcurve_ctrl) <- paste("Replicate_",unique(cfg_info$dataset_to_gradientpart$Replicate),sep="")
      colnames(norm_factors_meltingcurve_ctrl) <- cfg_info$temp_gradient
      norm_factors_meltingcurve_treat <- norm_factors_meltingcurve_ctrl
      
      if(ndata == 1){pdf("Meltingcurve normalization - with abundance effect.pdf")}
      if(ndata == 2){pdf("Meltingcurve normalization - without abundance effect.pdf")}
      
      if(class(fit) != "try-error")
      {
        plot(temps,relfcsmedian,main="Median melting curve",pch=20,xlab = "Temperature [?C]",ylab = "Fraction non-denatured")
        lines(x = new$x, y=predict(fit,newdata=new))
        
        ###relfcs at discrete temps of experiment on median curve fit
        relfc_curvefit_at_temps <- matrix(ncol=length(cfg_info$temp_gradient),nrow=1)
        for(t in cfg_info$temp_gradient)
        {
          ind <- which(cfg_info$temp_gradient == t)
          if(length(which(new$x == as.numeric(t))) > 0)###relfc at temperature can be directly taken from fit
          {
            relfc_curvefit_at_temps[1,ind] <- predict(fit,newdata=new)[which(new$x == as.numeric(t))]
          }else ###relfc as to be interpolated for temperature from fit
          {
            xlower <- which(new$x > as.numeric(t))[1]
            xhigher <- xlower+1
            relfc_curvefit_at_temps[1,ind] <- mean(predict(fit,newdata=new)[xlower:xhigher],na.rm=T)
          }
        }
        ###determine normalization factors for cond1 and cond2 for each temperature and each replicate
        for(r in unique(cfg_info$dataset_to_gradientpart$Replicate)) ###for each replicate
        {
          relfcsmedianctrl_currep <- colMedians(relfcsctrl[seq(r, length(relfcsctrl[,1]), length(unique(cfg_info$dataset_to_gradientpart$Replicate))),],na.rm=T)
          relfcsmediantreat_currep <- colMedians(relfcstreat[seq(r, length(relfcstreat[,1]), length(unique(cfg_info$dataset_to_gradientpart$Replicate))),],na.rm=T)
          norm_factors_meltingcurve_ctrl[r,] <- as.numeric(relfc_curvefit_at_temps)/relfcsmedianctrl_currep
          norm_factors_meltingcurve_treat[r,] <- as.numeric(relfc_curvefit_at_temps)/relfcsmediantreat_currep
          norm_factors_meltingcurve_ctrl[r,1] <- 1
          norm_factors_meltingcurve_treat[r,1] <- 1
          points(temps,relfcsmedianctrl_currep,pch=20,col=colors_melt_cond1[r])
          points(temps,relfcsmediantreat_currep,pch=20,col=colors_melt_cond2[r])
        }
        legend("topright", lty=1,lwd=2,col=c("black",colors_melt_cond1[unique(cfg_info$dataset_to_gradientpart$Replicate)],colors_melt_cond2[unique(cfg_info$dataset_to_gradientpart$Replicate)]), c("All",paste(cfg_info$condition_name_1,unique(cfg_info$dataset_to_gradientpart$Replicate),sep="_"),paste(cfg_info$condition_name_2,unique(cfg_info$dataset_to_gradientpart$Replicate),sep="_")), bty="o",cex=0.8, box.col="black",title = "Median per sample")
      }
      dev.off()
      if(ndata == 1)
      {
        norm_factors_meltingcurve_output[["norm_factors_meltingcurve_ctrl_with_abundance"]] <- norm_factors_meltingcurve_ctrl
        norm_factors_meltingcurve_output[["norm_factors_meltingcurve_treat_with_abundance"]] <- norm_factors_meltingcurve_treat
      }
      if(ndata == 2)
      {
        norm_factors_meltingcurve_output[["norm_factors_meltingcurve_ctrl_without_abundance"]] <- norm_factors_meltingcurve_ctrl
        norm_factors_meltingcurve_output[["norm_factors_meltingcurve_treat_without_abundance"]] <- norm_factors_meltingcurve_treat
      }
      norm_factors_meltingcurve <- rbind(norm_factors_meltingcurve_ctrl,norm_factors_meltingcurve_treat)
      ctrltreatcol <- as.data.frame(append(paste(cfg_info$condition_name_1,unique(cfg_info$dataset_to_gradientpart$Replicate),sep="_"),paste(cfg_info$condition_name_2,unique(cfg_info$dataset_to_gradientpart$Replicate),sep="_")))
      colnames(ctrltreatcol) <- "Condition"
      norm_factors_meltingcurve <- cbind(ctrltreatcol,norm_factors_meltingcurve)
      if(ndata == 1){saveExcel(Data = norm_factors_meltingcurve,File = "Meltingcurve normalization - with abundance effect.xlsx")}
      if(ndata == 2){saveExcel(Data = norm_factors_meltingcurve,File = "Meltingcurve normalization - without abundance effect.xlsx")}
    }
    
    ###now fit melting curves per protein with normalized relfcs
    
    if(progressbar == T){pb <- winProgressBar(title = paste("Fitting melting curves:",ndata),label=paste( round(0/(numreplicates*nrow(summary_data_combined[[ndata]]))*100, 0),"% done"), min = 0,max = numreplicates*nrow(summary_data_combined[[ndata]]), width = 300)}
    for(r in unique(cfg_info$dataset_to_gradientpart$Replicate)) ###for each replicate
    {
      for(i in 1:nrow(summary_data_combined[[ndata]])) ###each protein
      {
        ####summarize relfcs for cur protein and cur replicate for all temps
        relfcscond1 <- matrix(ncol=length(cfg_info$temp_gradient),nrow=1)
        relfcscond2 <- relfcscond1
        for(t in cfg_info$temp_gradient)
        {
          cols <- which(grepl(paste(t,"C_cond1.",r,sep=""),colnames(summary_data_combined[[ndata]])))
          ind <- which(cfg_info$temp_gradient == t)
          relfcscond1[1,ind] <- median(as.numeric(summary_data_combined[[ndata]][i,cols]),na.rm=T)
          cols <- which(grepl(paste(t,"C_cond2.",r,sep=""),colnames(summary_data_combined[[ndata]])))
          relfcscond2[1,ind] <- median(as.numeric(summary_data_combined[[ndata]][i,cols]),na.rm=T)
        }
        if(cfg_info$normalize_melting_curve == 1)
        {
          relfcscond1 <- as.numeric(relfcscond1)*norm_factors_meltingcurve_ctrl[r,]
          relfcscond2 <- as.numeric(relfcscond2)*norm_factors_meltingcurve_treat[r,]
        }else
        {
          relfcscond1 <- as.numeric(relfcscond1)
          relfcscond2 <- as.numeric(relfcscond2)
        }
        ####fit curves for cond1
        y <- relfcscond1
        y[1] <- median(as.numeric(summary_data_combined[[ndata]][i,which(grepl("37C_cond1.",colnames(summary_data_combined[[ndata]])))]),na.rm=T)
        fit <- fitSigmoidTR2(temps,y,startparam,10,y[1])
        new = data.frame(x = seq(min(temps),max(temps),len=500))
        if(class(fit) != "try-error")
        {
          if(ndata == 1)###data with abundance effect require a different melting point determination
          {
            inflec <- y[1] - ((y[1] - predict(fit,data.frame(x=100)))/2)
            try(Meltingpoints[i,r] <- (new$x[which(predict(fit,newdata=new)<=inflec)[1]]+new$x[which(predict(fit,newdata=new)<=inflec)[1]-1])/2,silent=TRUE)
          }else
          {
            pars <- coefficients(fit)
            a <- pars[["a"]]
            b <- pars[["b"]]
            pl <- pars[["Pl"]]
            Meltingpoints[i,r] <- a / (b - log((1-pl)/(1/2 - pl) - 1))
          }
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
          RSqares[i,r] <- (1-(sum(distres[1,]^2,na.rm=T)/sum(disttot[1,]^2,na.rm=T)))
        }
        ####fit curves for cond2
        ###check if there is a significant abundance change at reftemp, if yes take median of replicates as startpoint for all replicates, if not start at 1 as control
        # if(summary_data_combined[[ndata]]$pValue37C[i] < cfg_info$pvalue_abundance_cutoff & abs(round(2^summary_data_combined[[ndata]]$meanlog2Ratio37C[i]-1,2)) >= cfg_info$abundance_cutoff)
        # {
        #   y <- relfcscond2
        #   y[1] <- median(as.numeric(summary_data_combined[[ndata]][i,which(grepl("37C_cond2.",colnames(summary_data_combined[[ndata]])))]),na.rm=T)
        # }else ##no significant abundance change at reftemp -> start as Cond1 at 1
        # {
        #   y <- relfcscond2
        #   y[1] <- 1
        # }
        y <- relfcscond2
        y[1] <- median(as.numeric(summary_data_combined[[ndata]][i,which(grepl("37C_cond2.",colnames(summary_data_combined[[ndata]])))]),na.rm=T)
        fit <- fitSigmoidTR2(temps,y,startparam,10,y[1])
        new = data.frame(x = seq(min(temps),max(temps),len=500))
        if(class(fit) != "try-error")
        {
          if(ndata == 1)###data with abundance effect require a different melting point determination
          {
            inflec <- y[1] - ((y[1] - predict(fit,data.frame(x=100)))/2)
            try(Meltingpoints[i,numreplicates+r] <- (new$x[which(predict(fit,newdata=new)<=inflec)[1]]+new$x[which(predict(fit,newdata=new)<=inflec)[1]-1])/2,silent=TRUE)
          }else
          {
            pars <- coefficients(fit)
            a <- pars[["a"]]
            b <- pars[["b"]]
            pl <- pars[["Pl"]]
            Meltingpoints[i,numreplicates+r] <- a / (b - log((1-pl)/(1/2 - pl) - 1))
          }
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
          RSqares[i,numreplicates+r] <- (1-(sum(distres[1,]^2,na.rm=T)/sum(disttot[1,]^2,na.rm=T)))
        }
        counter <- counter + 1
        if(progressbar == T){setWinProgressBar(pb, counter, label=paste(round(counter/(numreplicates*nrow(summary_data_combined[[ndata]]))*100, 0)," % done (",counter,"/",numreplicates*nrow(summary_data_combined[[ndata]]),")",sep = ""))}
      }
    }
    if(progressbar == T){close(pb)}
    
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],Meltingpoints)
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],RSqares)
    
    #####add calculations for Melting curve Plateaus
    data <- summary_data_combined[[ndata]]
    plateaus <- as.data.frame(matrix(ncol=3,nrow=nrow(data)))
    colnames(plateaus) <- c("Plateau","Plateau_relfc_cond1","Plateau_relfc_cond2")
    max <- nrow(data)
    if(progressbar == T){pb <- winProgressBar(title = paste("Calculate plateaus:",ndata),label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)}
    Rsqcolscond1 <- which(grepl("Rsq_cond1",colnames(data)))
    Rsqcolscond2 <- which(grepl("Rsq_cond2",colnames(data)))
    for(i in 1:nrow(data))
    {
      if(length(which(!is.na(data[i,Rsqcolscond1]))) >= 1 & length(which(!is.na(data[i,Rsqcolscond2]))) >= 1)
      {
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
        ####calculate median at each temp
        relfcscond1 <- as.data.frame(t(colMedians(relfcscond1,na.rm=T)))
        relfcscond2 <- as.data.frame(t(colMedians(relfcscond2,na.rm=T)))
        if(cfg_info$normalize_melting_curve == 1)
        {
          ###normalize
          relfcscond1 <- relfcscond1[1:cfg_info$replicates,] * norm_factors_meltingcurve_ctrl
          relfcscond2 <- relfcscond2[1:cfg_info$replicates,] * norm_factors_meltingcurve_treat
        }
        x <- as.numeric(cfg_info$temp_gradient)
        ycond1 <- as.numeric(relfcscond1[1,])
        ycond2 <- as.numeric(relfcscond2[1,])
        fitcond1 <- fitSigmoidTR2(x,ycond1,startparam,10,ycond1[1])
        fitcond2 <- fitSigmoidTR2(x,ycond2,startparam,10,ycond2[1])
        new = data.frame(x = seq(min(x),max(x),len=500))
        if(class(fitcond1) != "try-error" & class(fitcond2) != "try-error") ###fitted curve in both
        {
          plateaucond1 <- predict(fitcond1,data.frame(x=100)) ###relfc plateau of melting curve of cond1s
          plateaucond2 <- predict(fitcond2,data.frame(x=100)) ###relfc plateau of melting curve of Cond2
          ####define temp from which on platea is reached in cond1 and cond2
          ###corresponds to temp where relfc < platea + 5 %
          ###use platea temp either of cond1 or cond2 depending which one is reaching platea later. use the later temp
          plateautempcond1 <- new$x[which(predict(fitcond1,newdata=new) < plateaucond1 + 0.1)[1]]
          plateautempcond2 <- new$x[which(predict(fitcond2,newdata=new) < plateaucond2 + 0.1)[1]]
          plateaus[i,1] <- max(plateautempcond1,plateautempcond2)
          plateaus[i,2] <- plateaucond1 + 0.1
          plateaus[i,3] <- plateaucond2 + 0.1
        }else ###no fit or one fit missing
        {
          plateaus[i,1] <- NA
          plateaus[i,2] <- NA
          plateaus[i,3] <- NA
        }
      }
      if(progressbar == T){setWinProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,")",sep = ""))}
    }
    if(progressbar == T){close(pb)}
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],plateaus)
    
    ####Mean delta Tm and significance of meltpoint difference (condition2 - condition1)
    Meltsignificance <- as.data.frame(matrix(ncol=2,nrow=nrow(summary_data_combined[[ndata]])))
    colnames(Meltsignificance) <- c("Mean_dTm","dTm_pValue")
    for(i in 1:nrow(summary_data_combined[[ndata]]))
    {
      Meltsignificance$Mean_dTm[i] <- mean(as.numeric(Meltingpoints[i,c((numreplicates+1):(2*numreplicates))]),na.rm=T) - mean(as.numeric(Meltingpoints[i,c(1:numreplicates)]),na.rm=T)
      Meltsignificance$dTm_pValue[i] <- ifelse(sum(!is.na(Meltingpoints[i,]))>=4,t.test2(as.numeric(Meltingpoints[i,c(1:numreplicates)]),as.numeric(Meltingpoints[i,c((numreplicates+1):(2*numreplicates))]), alternative = "two.sided", paired = FALSE,var.equal = TRUE),NA)
    }
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],Meltsignificance)
    summary_data_combined[[ndata]] <- summary_data_combined[[ndata]][order(summary_data_combined[[ndata]]$Protein),]
    summary_data_combined[[ndata]][is.na(summary_data_combined[[ndata]])] <- NA
  }
  output <- list()
  output[["summary_data_combined"]] <- summary_data_combined
  output[["norm_factors_meltingcurve"]] <- norm_factors_meltingcurve_output
  
  return(output)
}