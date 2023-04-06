t.test2 <- function(...) ###T-test
{
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

calculatepvalues = function(summary_data_combined,cfg_info,fordata = "both",progressbar=T) ###calculate pValues for log2 relfcs at each temperature
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
    p_valless <- as.data.frame(matrix(ncol=length(cfg_info$temp_gradient),nrow=nrow(summary_data_combined[[ndata]])))
    p_valgreater <- as.data.frame(matrix(ncol=length(cfg_info$temp_gradient),nrow=nrow(summary_data_combined[[ndata]])))
    max <- length(cfg_info$temp_gradient) * nrow(summary_data_combined[[ndata]])
    count <- 0
    if(progressbar == T){pb <- tcltk::tkProgressBar(title = paste("Calculating pValues:",ndata),label=paste( round(count/max*100, 0)," % done (",count,"/",max,")",sep = ""), min = 0,max = max, width = 300)}
    for(t in cfg_info$temp_gradient)
    {
      tindex <- which(cfg_info$temp_gradient == t)
      colnames(p_valless)[tindex] <- paste("pValue_less",t,"C",sep="")
      colnames(p_valgreater)[tindex] <- paste("pValue_greater",t,"C",sep="")
      cond1cols <- which(grepl(paste(t,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))
      cond2cols <- which(grepl(paste(t,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))
      for(i in 1:nrow(summary_data_combined[[ndata]]))
      {
        ####remove relfc values for pvalue calculation if both cond2 and cond1 is < cutoff
        tempcond1s <- summary_data_combined[[ndata]][i,cond1cols]
        tempcond2s <- summary_data_combined[[ndata]][i,cond2cols]
        
        if(cfg_info$remove_weak_quant == 1)
        {
          for(q in 1:length(cond1cols))
          {
            if(!is.na(tempcond1s[q]) & !is.na(tempcond2s[q]))
            {
              if(tempcond1s[q] < cfg_info$relfc_cutoff_plot & tempcond2s[q] < cfg_info$relfc_cutoff_plot)
              {
                tempcond1s[q] <- as.numeric(NA)
                tempcond2s[q] <- as.numeric(NA)
              }
            }else
            {
              tempcond1s[q] <- as.numeric(NA)
              tempcond2s[q] <- as.numeric(NA)
            }
          }
        }
        
        ####need at least 2 replicate values for cond1 and cond2 so in total 4 values
        if(tindex == 1)###one sample t-test @ reference temperature
        {
          p_valgreater[i,tindex] = ifelse(sum(!is.na(tempcond2s))>=2,t.test2(log2(tempcond2s),mu=0, alternative = "greater", paired = FALSE,var.equal = TRUE),NA)
          p_valless[i,tindex] = ifelse(sum(!is.na(tempcond2s))>=2,t.test2(log2(tempcond2s),mu=0, alternative = "less", paired = FALSE,var.equal = TRUE),NA)
        }else
        {
          tempcond2s <- tempcond2s/tempcond1s
          p_valgreater[i,tindex] = ifelse(sum(!is.na(tempcond2s))>=2,t.test2(log2(tempcond2s),mu=0, alternative = "greater", paired = FALSE,var.equal = TRUE),NA)
          p_valless[i,tindex] = ifelse(sum(!is.na(tempcond2s))>=2,t.test2(log2(tempcond2s),mu=0, alternative = "less", paired = FALSE,var.equal = TRUE),NA)
        }
        
        count <- count + 1
        if(progressbar == T){tcltk::setTkProgressBar(pb, count, label=paste( round(count/max*100, 0)," % done (",count,"/",max,")",sep = ""))}
      }
    }
    if(progressbar == T){close(pb)}
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],p_valless)
    summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],p_valgreater)
  }
  
  ####plot density distribution of pvalues
  
  for(ndata in calcsetup)
  {
    pdf(paste("pvalue distributions ",ndata,".pdf",sep=""))
    for(temp in cfg_info$temp_gradient) ###for each temperature
    {
      if(ndata == 2 & temp != "37" | ndata == 1)
      {
        pvals <- subset(summary_data_combined[[ndata]], select = c(paste("pValue_less",temp,"C",sep=""),paste("pValue_greater",temp,"C",sep="")))
        if(any(!is.na(pvals[,1])))plot(density(as.matrix(subset(summary_data_combined[[ndata]], select = c(paste("pValue_less",temp,"C",sep="")))),na.rm=T),main=paste("pvalues_less at ",temp,"°C",sep=""),xlab=paste("pvalue ",temp,"°C",sep=""),col="black")
        if(any(!is.na(pvals[,2])))plot(density(as.matrix(subset(summary_data_combined[[ndata]], select = c(paste("pValue_greater",temp,"C",sep="")))),na.rm=T),main=paste("pvalues_greater at ",temp,"°C",sep=""),xlab=paste("pvalue ",temp,"°C",sep=""),col="black")                    
      }
    }
    dev.off()
  }
  
  return(summary_data_combined)
}