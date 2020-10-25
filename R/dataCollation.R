dataCollation = function(summary_data_combined,cfg_info,fordata = "both",progressbar = T) ###collate log2ratios and pvalues
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
    if(cfg_info$pVal_collation_function == 3) ### collate by PECA
    {
      #######calculate significant summed log2 ratios and p-Values#########
      collationlog2ratio <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValue <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationcount <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationstart <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      #prepare matrix on which PECA will be applied
      peca_matrix <- as.data.frame(matrix(ncol=1+(cfg_info$replicates*2),nrow=nrow(summary_data_combined[[ndata]])*length(cfg_info$temp_gradient)))
      colnames(peca_matrix) <- c(paste("Cond1_",1:cfg_info$replicates,sep=""),paste("Cond2_",1:cfg_info$replicates,sep=""),"Protein")
      peca_matrix_index_counter <- 0
      for(c in 1:(ncol(peca_matrix)-1))
      {
        peca_matrix[,c] <- as.numeric(peca_matrix[,c])
      }
      peca_matrix[,ncol(peca_matrix)] <- as.character(peca_matrix[,ncol(peca_matrix)])
      
      if(progressbar == T){pb <- winProgressBar(title = paste("Collating data:",ndata),label=paste( round(0/nrow(summary_data_combined[[ndata]])*100, 0),"% done"), min = 0,max = nrow(summary_data_combined[[ndata]]), width = 300)}
      for(i in 1:nrow(summary_data_combined[[ndata]]))
      {
        sum <- 0
        pValue <- NA
        finalpValue <- NA
        collationcounter <- 0
        counter <- 0
        direction <- 0 ### -1 if lower ratio < 0 and 1 if ratio > 0
        relevanttemps <- NA
        Values <- matrix(nrow = 0,ncol = 3) ### groups for data collation are summarized here to later find the optimum
        colnames(Values) <- c("log2diff","count","start")
        
        temp_peca_matrix <- as.data.frame(matrix(ncol=2+(cfg_info$replicates*2),nrow=length(cfg_info$temp_gradient)))
        temp_peca_matrix_index_counter <- 0
        for(c in 1:(ncol(temp_peca_matrix)-1))
        {
          temp_peca_matrix[,c] <- as.numeric(temp_peca_matrix[,c])
        }
        temp_peca_matrix[,ncol(temp_peca_matrix)] <- as.character(temp_peca_matrix[,ncol(peca_matrix)])
        temp_peca_matrix_index_counter <- 0
        temp_peca_matrix_block_counter <- 0
        
        for(temp in cfg_info$temp_gradient) 
        {
          test2 <- subset(summary_data_combined[[ndata]], select = c(paste("meanlog2Ratio",temp,"C",sep=""),paste("pValue_less",temp,"C",sep=""),paste("pValue_greater",temp,"C",sep="")))[i,]
          minindex <- 1+which(test2[1,2:3] == min(test2[1,2:3],na.rm=T))
          minindex <- ifelse(length(minindex)>0,minindex,2)
          test <- t(as.data.frame(c(test2[1,1],test2[1,minindex])))
          colnames(test) <- c("log2diff","pValue")
          rownames(test) <- c()
          
          if(!is.na(test[1]) && !is.na(test[2]) && test[2] <= cfg_info$pVal_collation_cutoff)
          {
            if(counter == 0) #Zähle wieviele Werte in einer Reihe zusammen unter dem Cutoff liegen
            {
              if(test[1] > 0){direction <- 1}
              if(test[1] < 0){direction <- -1}     
              collationcounter <- collationcounter + 1
              counter <- 1
              pValue <- test[2]
              relevanttemps <- temp
              sum <- test[1]
            }else
            {
              if(direction == 1 && test[1]>0 | direction == -1 && test[1]<0)
              {
                counter <- counter + 1
                pValue <- rbind(pValue,test[2])
                relevanttemps <- rbind(relevanttemps,temp)
                sum <- sum + test[1]
              }else ###significant difference but not in the correct direction then stop the combination here and start a new one from this point
              {
                ###combine statistics using PECA
                #prepare raw unlogged data for aggregation
                #samples in cols, temperatures in rows
                datamatrix <- matrix(nrow = 0,ncol = 2*cfg_info$replicates)
                for(relevanttemp in relevanttemps)
                {
                  cond1cols <- which(grepl(paste(relevanttemp,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                  cond2cols <- which(grepl(paste(relevanttemp,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                  datamatrix <- rbind(datamatrix,as.numeric((summary_data_combined[[ndata]])[i,c(cond1cols,cond2cols)]))
                }
                datamatrix <- as.data.frame(datamatrix*100)
                temp_peca_matrix_block_counter <- temp_peca_matrix_block_counter + 1
                datamatrix$block <- temp_peca_matrix_block_counter
                datamatrix$id <- summary_data_combined[[ndata]]$Protein[i]
                
                #write into global matrix
                start <- temp_peca_matrix_index_counter+1
                end <- temp_peca_matrix_index_counter + nrow(datamatrix)
                set(temp_peca_matrix,i = as.integer(start:end),j=as.integer(1:((2*cfg_info$replicates)+1)),datamatrix[,1:((2*cfg_info$replicates)+1)])
                set(temp_peca_matrix,i = as.integer(start:end),j=as.integer(ncol(temp_peca_matrix)),datamatrix[,ncol(datamatrix)])
                temp_peca_matrix_index_counter <- temp_peca_matrix_index_counter + nrow(datamatrix)

                
                Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
                colnames(Values) <- c("log2diff","count","start")
                ######now go a step back and do a new collation at this point
                if(test[1] > 0){direction <- 1}
                if(test[1] < 0){direction <- -1}     
                collationcounter <- collationcounter + 1
                counter <- 1
                pValue <- test[2]
                relevanttemps <- temp
                sum <- test[1]
              }
            }  
          }else if(counter != 0) ###if there was one value below the cutoff then save the pvalue and sumlog2diff of that collation and reset counter for next round of collation
          {
            ###combine statistics using PECA
            #prepare raw unlogged data for aggregation
            #samples in cols, temperatures in rows
            datamatrix <- matrix(nrow = 0,ncol = 2*cfg_info$replicates)
            for(relevanttemp in relevanttemps)
            {
              cond1cols <- which(grepl(paste(relevanttemp,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
              cond2cols <- which(grepl(paste(relevanttemp,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
              datamatrix <- rbind(datamatrix,as.numeric((summary_data_combined[[ndata]])[i,c(cond1cols,cond2cols)]))
            }
            datamatrix <- as.data.frame(datamatrix*100)
            temp_peca_matrix_block_counter <- temp_peca_matrix_block_counter + 1
            datamatrix$block <- temp_peca_matrix_block_counter
            datamatrix$id <- summary_data_combined[[ndata]]$Protein[i]
            
            #write into global matrix
            start <- temp_peca_matrix_index_counter+1
            end <- temp_peca_matrix_index_counter + nrow(datamatrix)
            set(temp_peca_matrix,i = as.integer(start:end),j=as.integer(1:((2*cfg_info$replicates)+1)),datamatrix[,1:((2*cfg_info$replicates)+1)])
            set(temp_peca_matrix,i = as.integer(start:end),j=as.integer(ncol(temp_peca_matrix)),datamatrix[,ncol(datamatrix)])
            temp_peca_matrix_index_counter <- temp_peca_matrix_index_counter + nrow(datamatrix)
            
            Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
            colnames(Values) <- c("log2diff","count","start")
            counter <- 0
          }
          if(counter != 0 && which(cfg_info$temp_gradient == temp) == length(cfg_info$temp_gradient)) ###if also for the last temperature a significant difference was observed
          {
            ###combine statistics using PECA
            #prepare raw unlogged data for aggregation
            #samples in cols, temperatures in rows
            datamatrix <- matrix(nrow = 0,ncol = 2*cfg_info$replicates)
            for(relevanttemp in relevanttemps)
            {
              cond1cols <- which(grepl(paste(relevanttemp,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
              cond2cols <- which(grepl(paste(relevanttemp,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
              datamatrix <- rbind(datamatrix,as.numeric((summary_data_combined[[ndata]])[i,c(cond1cols,cond2cols)]))
            }
            datamatrix <- as.data.frame(datamatrix*100)
            temp_peca_matrix_block_counter <- temp_peca_matrix_block_counter + 1
            datamatrix$block <- temp_peca_matrix_block_counter
            datamatrix$id <- summary_data_combined[[ndata]]$Protein[i]
            
            #write into global matrix
            start <- temp_peca_matrix_index_counter+1
            end <- temp_peca_matrix_index_counter + nrow(datamatrix)
            set(temp_peca_matrix,i = as.integer(start:end),j=as.integer(1:((2*cfg_info$replicates)+1)),datamatrix[,1:((2*cfg_info$replicates)+1)])
            set(temp_peca_matrix,i = as.integer(start:end),j=as.integer(ncol(temp_peca_matrix)),datamatrix[,ncol(datamatrix)])
            temp_peca_matrix_index_counter <- temp_peca_matrix_index_counter + nrow(datamatrix)
            
            
            Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
            colnames(Values) <- c("log2diff","count","start")
            counter <- 0
          }
        }
        Values <- as.data.frame(Values)
        
        Values$log2diff <- as.numeric(as.character(Values$log2diff))
        Values$count <- as.numeric(as.character(Values$count))
        Values$start <- as.numeric(as.character(Values$start))
        Values <- subset(Values,!is.na(log2diff))
        
        ##### select Group with abs(max log2ratio)
        selectedsum <- NA
        selectedpVal <- NA
        selectedcount <- NA
        selectedstart <- NA
        if(nrow(Values)>0)
        {
          sel <- which(abs(Values$log2diff) == max(abs(Values$log2diff)))
          
          selectedsum <- Values$log2diff[sel]
          selectedcount <- Values$count[sel]
          selectedstart <- Values$start[sel]
          
          ##store respective matrix for PECA
          datamatrix <- temp_peca_matrix[which(temp_peca_matrix[,ncol(temp_peca_matrix)-1] == as.numeric(rownames(Values)[sel])),-(ncol(temp_peca_matrix)-1)]
          start <- peca_matrix_index_counter+1
          end <- peca_matrix_index_counter + nrow(datamatrix)
          set(peca_matrix,i = as.integer(start:end),j=as.integer(1:((2*cfg_info$replicates))),datamatrix[,1:((2*cfg_info$replicates))])
          set(peca_matrix,i = as.integer(start:end),j=as.integer(ncol(peca_matrix)),datamatrix[,ncol(peca_matrix)])
          peca_matrix_index_counter <- peca_matrix_index_counter + nrow(datamatrix)

        }else ### if the pvalue cutoff was never reached, then use the values at abs(maxlog2diff)
        {
          ####find temp with max meanlog2Ratio
          tmp <- summary_data_combined[[ndata]][i,which(grepl("meanlog2Ratio",colnames(summary_data_combined[[ndata]])))]
          tmptemp <- cfg_info$temp_gradient[which(abs(tmp) == max(abs(tmp),na.rm=T))] ###temp at max ratio
          if(length(tmptemp)>0)
          {
            selectedsum <- summary_data_combined[[ndata]][i,paste("meanlog2Ratio",tmptemp[1],"C",sep="")]
            selectedcount <- 1
            selectedstart <- which(cfg_info$temp_gradient == tmptemp[1])
            
            datamatrix <- matrix(nrow = 0,ncol = 2*cfg_info$replicates)
            for(relevanttemp in tmptemp)
            {
              cond1cols <- which(grepl(paste(relevanttemp,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
              cond2cols <- which(grepl(paste(relevanttemp,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
              datamatrix <- rbind(datamatrix,as.numeric((summary_data_combined[[ndata]])[i,c(cond1cols,cond2cols)]))
            }
            datamatrix <- as.data.frame(datamatrix*100)
            temp_peca_matrix_block_counter <- temp_peca_matrix_block_counter + 1
            datamatrix$id <- summary_data_combined[[ndata]]$Protein[i]
            
            start <- peca_matrix_index_counter+1
            end <- peca_matrix_index_counter + nrow(datamatrix)
            set(peca_matrix,i = as.integer(start:end),j=as.integer(1:((2*cfg_info$replicates))),datamatrix[,1:((2*cfg_info$replicates))])
            set(peca_matrix,i = as.integer(start:end),j=as.integer(ncol(peca_matrix)),datamatrix[,ncol(peca_matrix)])
            peca_matrix_index_counter <- peca_matrix_index_counter + nrow(datamatrix)
          }
          
        }
        
        collationlog2ratio[i,1] <- selectedsum
        collationcount[i,1] <- selectedcount
        collationstart[i,1] <- selectedstart
        
        if(progressbar == T){setWinProgressBar(pb, i, label=paste( round(i/nrow(summary_data_combined[[ndata]])*100, 0)," % done (",i,"/",nrow(summary_data_combined[[ndata]]),")",sep = ""))}
      }
      if(progressbar == T){close(pb)}
      
      ##Now perform PECA
      peca_res <- PECA_df(peca_matrix,id = "Protein",samplenames1 = paste("Cond2_",1:cfg_info$replicates,sep=""),samplenames2 = paste("Cond1_",1:cfg_info$replicates,sep=""),test = "t")
      peca_res <- peca_res[match(summary_data_combined[[ndata]]$Protein,rownames(peca_res)),]
      
      collationlog2ratio[is.infinite(collationlog2ratio$V1),1] <- NA
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationlog2ratio)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.log2Ratio"
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],peca_res$p)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue"
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationcount)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.count"
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationstart)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.start"
    }
    
    if(cfg_info$pVal_collation_function == 2) ### collate all available pvalues
    {
      collationlog2ratio <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      
      collationpValuelog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValuelesslog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValuegreaterlog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpvaluedirectionlog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpvaluecountlog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      
      collationratiocount <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationratiostart <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      if(progressbar == T){pb <- winProgressBar(title = paste("Collating data:",ndata),label=paste( round(0/nrow(summary_data_combined[[ndata]])*100, 0),"% done"), min = 0,max = nrow(summary_data_combined[[ndata]]), width = 300)}
      for(i in 1:nrow(summary_data_combined[[ndata]]))###for each protein
      {
        ###combine pvalues for less and greater, then decide which direction is most significant
        combinedpvaluelog2 <- NULL
        pvaluescountlog2 <- NULL
        
        for(d in c("less","greater"))
        {
          pvalues <- NULL
          rawdatalog2 <- NULL
          for(temp in cfg_info$temp_gradient) ###collect all pvalues and raw data for subsequent combination
          {
            test <- subset(summary_data_combined[[ndata]], select = c(paste("meanlog2Ratio",temp,"C",sep=""),paste("pValue_",d,temp,"C",sep="")))[i,]
            colnames(test) <- c("log2diff","pValue")
            testdatalog2 <- summary_data_combined[[ndata]][i,grepl(paste("^log2Ratio",temp,"C",sep=""),colnames(summary_data_combined[[ndata]]))]
            if(!is.na(test[1]) && !is.na(test[2]))
            {
              pvalues <- append(pvalues,as.numeric(test[2]))
              
              temps <- as.matrix(t(testdatalog2))
              ###more replicates for 37 °C so use a subset of the data as raw data input
              if(temp == "37")
              {
                selection <- c(1:cfg_info$replicates)
                temps <- temps[selection]
              }
              rownames(temps) <- c()
              rawdatalog2 <- cbind(rawdatalog2,temps)
            }
          }
          ###now combine pvalues
          if(length(pvalues) > 1) ###more than one pvalue
          {
            ###remove rows with missing values log2
            if(length(which(rowSums(is.na(rawdatalog2))==0)) == 1)
            {
              rawdatalog2 <- t(as.data.frame(rawdatalog2[rowSums(is.na(rawdatalog2))==0,]))
            }else if(length(which(rowSums(is.na(rawdatalog2))==0)) > 1)
            {
              rawdatalog2 <- as.data.frame(rawdatalog2[rowSums(is.na(rawdatalog2))==0,])
            }else
            {
              rawdatalog2 <- matrix(ncol=1,nrow=0)
            }
            
            
            if(nrow(rawdatalog2) > 1)
            {
              ###rawdata format --> nrow = length(pvalues) ncol = replicates or sample size
              combinedpvaluelog2 <- append(combinedpvaluelog2,empiricalBrownsMethod(data_matrix=t(rawdatalog2), p_values=pvalues, extra_info=FALSE))
              pvaluescountlog2 <- append(pvaluescountlog2,length(pvalues))
            }else###if not enough data to infere covariance matrix simply use the most significant pvalue
            {
              combinedpvaluelog2 <- append(combinedpvaluelog2,min(pvalues,na.rm=T))
              pvaluescountlog2 <- append(pvaluescountlog2,1)
            }
            
            
          }else if(length(pvalues) == 1) ###only one pvalue
          {
            combinedpvaluelog2 <- append(combinedpvaluelog2,min(pvalues,na.rm=T))
            pvaluescountlog2 <- append(pvaluescountlog2,1)
          }else
          {
            combinedpvaluelog2 <- append(combinedpvaluelog2,NA)
            pvaluescountlog2 <- append(pvaluescountlog2,0)
          }
        }
        if(length(which(!is.na(combinedpvaluelog2))) > 0)
        {
          collationpValuelesslog2[i,1] <- combinedpvaluelog2[1]
          collationpValuegreaterlog2[i,1] <- combinedpvaluelog2[2]
          collationpvaluedirectionlog2[i,1] <- c("less","greater")[which(combinedpvaluelog2 == min(combinedpvaluelog2,na.rm=T))[1]]
          collationpvaluecountlog2[i,1] <- pvaluescountlog2[which(combinedpvaluelog2 == min(combinedpvaluelog2,na.rm=T))[1]]
          
          
          ####now combine consecutive significant ratios with same direction as indicated by selected combined pvalues 
          sum <- 0
          collationcounter <- 0
          counter <- 0
          direction <- 0 ### -1 if lower ratio < 0 and 1 if ratio > 0
          relevanttemps <- NA
          Values <- matrix(nrow = 0,ncol = 3) ### groups for data collation are summarized here to later find the optimum
          colnames(Values) <- c("log2diff","count","start")
          for(temp in cfg_info$temp_gradient) 
          {
            test <- subset(summary_data_combined[[ndata]], select = c(paste("meanlog2Ratio",temp,"C",sep=""),paste("pValue_less",temp,"C",sep=""),paste("pValue_greater",temp,"C",sep="")))[i,]
            colnames(test) <- c("log2diff","pValue_less","pValue_greater")
            
            if(!is.na(test[1]) && !is.na(test[2]) && test[2] <= cfg_info$pVal_collation_cutoff | !is.na(test[1]) && !is.na(test[3]) && test[3] <= cfg_info$pVal_collation_cutoff)
            {
              if(counter == 0) #Zähle wieviele Werte in einer Reihe zusammen unter dem Cutoff liegen
              {
                if(test[1] > 0){direction <- 1}
                if(test[1] < 0){direction <- -1}     
                collationcounter <- collationcounter + 1
                counter <- 1
                relevanttemps <- temp
                sum <- test[1]
              }else
              {
                if(direction == 1 && test[1]>0 | direction == -1 && test[1]<0)
                {
                  counter <- counter + 1
                  relevanttemps <- rbind(relevanttemps,temp)
                  sum <- sum + test[1]
                }else ###significant difference but not in the correct direction then stop the combination here and start a new one from this point
                {
                  Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
                  colnames(Values) <- c("log2diff","count","start")
                  ######now go a step back and do a new collation at this point
                  if(test[1] > 0){direction <- 1}
                  if(test[1] < 0){direction <- -1}     
                  collationcounter <- collationcounter + 1
                  counter <- 1
                  relevanttemps <- temp
                  sum <- test[1]
                }
              }  
            }else if(counter != 0) ###if there was one value below the cutoff then save the pvalue and sumlog2diff of that collation and reset counter for next round of collation
            {
              Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
              colnames(Values) <- c("log2diff","count","start")
              counter <- 0
            }
            if(counter != 0 && which(cfg_info$temp_gradient == temp) == length(cfg_info$temp_gradient)) ###if also for the last temperature a significant difference was observed
            {
              Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
              colnames(Values) <- c("log2diff","count","start")
              counter <- 0
            }
          }
          Values <- as.data.frame(Values)
          Values <- subset(Values,!is.na(log2diff))
          ##### select Group with abs(max log2ratio)
          selectedsum <- NA
          selectedcount <- NA
          selectedstart <- NA
          if(nrow(Values)>0)
          {
            max <- Values[which.max(Values$log2diff),]
            min <- Values[which.min(Values$log2diff),]
            if(abs(unlist(max$log2diff)) >= abs(unlist(min$log2diff)))
            {
              selectedsum <- unlist(max$log2diff)
              selectedcount <- unlist(max$count)
              selectedstart <- unlist(max$start)
            }else
            {
              selectedsum <- unlist(min$log2diff)
              selectedcount <- unlist(min$count)
              selectedstart <- unlist(min$start)
            } 
          }else ### if the pvalue cutoff was never reached, then use the value at most significant temp in the direction of selected collated pvalues
          {
            testdata <- summary_data_combined[[ndata]][i,which(grepl("meanlog2Ratio",colnames(summary_data_combined[[ndata]])))]
            testpvalues <- summary_data_combined[[ndata]][i,which(grepl(paste("pValue_",collationpvaluedirectionlog2[i,1],sep=""),colnames(summary_data_combined[[ndata]])))]
            selectedsum <- testdata[1,which(testpvalues == min(testpvalues,na.rm=T))]
            selectedcount <- 1
            selectedstart <- which(testpvalues == min(testpvalues,na.rm=T)) ###temp at max ratio
          }
          
          collationlog2ratio[i,1] <- selectedsum
          collationratiocount[i,1] <- selectedcount
          collationratiostart[i,1] <- selectedstart
          
          if(!is.na(selectedsum))
          {
            if(selectedsum > 0)
            {
              collationpValuelog2[i,1] <- combinedpvaluelog2[2]
            }else
            {
              collationpValuelog2[i,1] <- combinedpvaluelog2[1]
            }
          }
          
        }
        if(progressbar == T){setWinProgressBar(pb, i, label=paste( round(i/nrow(summary_data_combined[[ndata]])*100, 0)," % done (",i,"/",nrow(summary_data_combined[[ndata]]),")",sep = ""))}
      }
      if(progressbar == T){close(pb)} 
      
      pdf(paste("Collated pvalue",ndata,".pdf",sep=""))
      plot(density(as.matrix(collationpValuelesslog2),na.rm=T),main="Collation pvalue_less")
      hist(as.numeric(as.matrix(collationpValuelesslog2)))
      plot(density(as.matrix(collationpValuegreaterlog2),na.rm=T),main="Collation pvalue_greater")
      hist(as.numeric(as.matrix(collationpValuegreaterlog2)))
      dev.off()
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationlog2ratio)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.log2Ratio"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpValuelog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpValuelesslog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue.less"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpValuegreaterlog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue.greater"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpvaluedirectionlog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue.direction"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpvaluecountlog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pvalue.count"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationratiocount)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.count"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationratiostart)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.start"
      
    }
    
    if(cfg_info$pVal_collation_function == 1) ### collate pvalues with a sliding-window approach
    {
      collationlog2ratio <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      
      collationpValuelog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValuelesslog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValuegreaterlog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValuelessstart <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValuegreaterstart <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpvaluedirectionlog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpvaluecountlog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpvaluestartlog2 <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      
      collationratiocount <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationratiostart <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      nGradientSlices <- 3 ###defines how many pvalues are collated per slice test e.g 9 temperatures in 3 slices --> 3 pvalues per slice are tested
      npvalsperslice <- floor(length(cfg_info$temp_gradient)/nGradientSlices)
      
      if(progressbar == T){pb <- winProgressBar(title = paste("Collating data:",ndata),label=paste( round(0/nrow(summary_data_combined[[ndata]])*100, 0),"% done"), min = 0,max = nrow(summary_data_combined[[ndata]]), width = 300)}
      for(i in 1:nrow(summary_data_combined[[ndata]]))###for each protein
      {
        ###combine pvalues for less and greater, then decide which direction is most significant
        combinedpvaluelog2 <- NULL
        pvaluescountlog2 <- NULL
        combinedpvaluestart <- NULL
        
        for(d in c("less","greater"))
        {
          pvalues <- NULL
          rawdatalog2 <- NULL
          for(temp in cfg_info$temp_gradient) ###collect all pvalues and raw data for subsequent combination
          {
            test <- subset(summary_data_combined[[ndata]], select = c(paste("meanlog2Ratio",temp,"C",sep=""),paste("pValue_",d,temp,"C",sep="")))[i,]
            colnames(test) <- c("log2diff","pValue")
            cond1relfcs <- summary_data_combined[[ndata]][i,grepl(paste("^",temp,"C_cond1",sep=""),colnames(summary_data_combined[[ndata]]))]
            cond2relfcs <- summary_data_combined[[ndata]][i,grepl(paste("^",temp,"C_cond2",sep=""),colnames(summary_data_combined[[ndata]]))]
            if(any(cond1relfcs >= cfg_info$relfc_cutoff_plot,na.rm = T) | any(cond2relfcs >= cfg_info$relfc_cutoff_plot,na.rm = T))
            {
              testdatalog2 <- log2(cond2relfcs/cond1relfcs)
            }else
            {
              testdatalog2 <- matrix(ncol=cfg_info$replicates,nrow=1)
            }
            if(!is.na(test[1]) && !is.na(test[2]))
            {
              pvalues <- append(pvalues,as.numeric(test[2]))
              
              temps <- as.matrix(t(testdatalog2))
              ###more replicates for 37 °C so use a subset of the data as raw data input
              if(temp == "37")
              {
                selection <- c(1:cfg_info$replicates)
                temps <- temps[selection]
              }
              rownames(temps) <- c()
              rawdatalog2 <- cbind(rawdatalog2,temps)
            }else
            {
              pvalues <- append(pvalues,NA)
              temps <- matrix(ncol=1,nrow=cfg_info$replicates)
              rownames(temps) <- c()
              rawdatalog2 <- cbind(rawdatalog2,temps)
            }
          }
          ###now combine pvalues
          if(length(!is.na(pvalues)) > 2) ###more than one pvalue
          {
            collatedpvalueslices <- NULL
            collatedpvalueslicescount <- NULL
            collatedpvalueslicesstart <- NULL
            ####Slice pvalue in slices containing nTemps/nslices pvalues
            for(s in 1:(length(cfg_info$temp_gradient)-npvalsperslice+1))
            {
              temppvalues <- pvalues[s:(s+npvalsperslice-1)]
              temprawdatalog2sliced <- rawdatalog2[,s:(s+npvalsperslice-1)]
              
              if(any(is.na(temppvalues)))
              {
                temprawdatalog2sliced <- temprawdatalog2sliced[,-c(which(is.na(temppvalues)))]
                temppvalues <- temppvalues[-c(which(is.na(temppvalues)))]
              }
              if(length(temppvalues) >= 2)
              {
                temprawdatalog2slicedsave <- temprawdatalog2sliced
                temppvaluessave <- temppvalues
                ###remove rows with missing values
                if(length(which(rowSums(is.na(temprawdatalog2sliced))==0)) == 1)
                {
                  temprawdatalog2sliced <- t(as.data.frame(temprawdatalog2sliced[rowSums(is.na(temprawdatalog2sliced))==0,]))
                }else if(length(which(rowSums(is.na(temprawdatalog2sliced))==0)) > 1)
                {
                  temprawdatalog2sliced <- as.data.frame(temprawdatalog2sliced[rowSums(is.na(temprawdatalog2sliced))==0,])
                }else
                {
                  temprawdatalog2sliced <- matrix(ncol=1,nrow=0)
                }
                
                if(nrow(temprawdatalog2sliced) > 1)
                {
                  ###rawdata format --> nrow = length(pvalues) ncol = replicates or sample size
                  collatedpvalueslices <- append(collatedpvalueslices,empiricalBrownsMethod(data_matrix=t(temprawdatalog2sliced), p_values=temppvalues, extra_info=FALSE))
                  collatedpvalueslicescount <- append(collatedpvalueslicescount,length(temppvalues))
                  collatedpvalueslicesstart <- append(collatedpvalueslicesstart,s)
                }else###if not enough data to infere covariance matrix 
                {
                  ###check if only one temperature contained NAs and in that case remove the pvalue and the raw data and check again
                  temprawdatalog2sliced <- temprawdatalog2slicedsave
                  
                  if(length(which(colSums(is.na(temprawdatalog2sliced))>0))>0)
                  {
                    temppvalues <- temppvalues[-which(colSums(is.na(temprawdatalog2sliced))>0)]
                    temprawdatalog2sliced <- temprawdatalog2sliced[,-which(colSums(is.na(temprawdatalog2sliced))>0)]
                  }
                  if(length(temppvalues) > 1)
                  {
                    if(nrow(temprawdatalog2sliced) > 1)
                    {
                      collatedpvalueslices <- append(collatedpvalueslices,empiricalBrownsMethod(data_matrix=t(temprawdatalog2sliced), p_values=temppvalues, extra_info=FALSE))
                      collatedpvalueslicescount <- append(collatedpvalueslicescount,length(temppvalues))
                      collatedpvalueslicesstart <- append(collatedpvalueslicesstart,s)
                    }else
                    {
                      temppvalues <- temppvaluessave
                      ###otherwise simply use the most significant pvalue
                      if(any(!is.na(temppvalues)))
                      {
                        collatedpvalueslices <- append(collatedpvalueslices,min(temppvalues,na.rm=T))
                        collatedpvalueslicescount <- append(collatedpvalueslicescount,1)
                        collatedpvalueslicesstart <- append(collatedpvalueslicesstart,(s-1+which(temppvalues==min(temppvalues,na.rm=T))))
                      }
                    }
                  }else
                  {
                    temppvalues <- temppvaluessave
                    ###otherwise simply use the most significant pvalue
                    if(any(!is.na(temppvalues)))
                    {
                      collatedpvalueslices <- append(collatedpvalueslices,min(temppvalues,na.rm=T))
                      collatedpvalueslicescount <- append(collatedpvalueslicescount,1)
                      collatedpvalueslicesstart <- append(collatedpvalueslicesstart,(s-1+which(temppvalues==min(temppvalues,na.rm=T))))
                    }
                  }
                  
                }
              }else
              {
                if(any(!is.na(temppvalues)))
                {
                  collatedpvalueslices <- append(collatedpvalueslices,min(temppvalues,na.rm=T))
                  collatedpvalueslicescount <- append(collatedpvalueslicescount,1)
                  collatedpvalueslicesstart <- append(collatedpvalueslicesstart,(s-1+which(temppvalues==min(temppvalues,na.rm=T))))
                }
              }
              
            }
            
            if(any(!is.na(collatedpvalueslices)))
            {
              combinedpvaluelog2 <- append(combinedpvaluelog2,min(collatedpvalueslices,na.rm=T))
              pvaluescountlog2 <- append(pvaluescountlog2,collatedpvalueslicescount[which(collatedpvalueslices == min(collatedpvalueslices,na.rm=T))[1]])
              combinedpvaluestart <- append(combinedpvaluestart,collatedpvalueslicesstart[which(collatedpvalueslices == min(collatedpvalueslices,na.rm=T))[1]])
            }
          }else if(length(pvalues) >= 1) ###only one pvalue
          {
            combinedpvaluelog2 <- append(combinedpvaluelog2,min(pvalues,na.rm=T))
            pvaluescountlog2 <- append(pvaluescountlog2,1)
            combinedpvaluestart <- append(pvaluescountlog2,which(pvalues==min(pvalues,na.rm=T))[1])
          }else
          {
            combinedpvaluelog2 <- append(combinedpvaluelog2,NA)
            pvaluescountlog2 <- append(pvaluescountlog2,0)
            combinedpvaluestart <- NA
          }
        }
        if(length(which(!is.na(combinedpvaluelog2))) > 0)
        {
          collationpValuelesslog2[i,1] <- combinedpvaluelog2[1]
          collationpValuegreaterlog2[i,1] <- combinedpvaluelog2[2]
          collationpValuelessstart[i,1] <- combinedpvaluestart[1]
          collationpValuegreaterstart[i,1] <- combinedpvaluestart[2]
          collationpvaluedirectionlog2[i,1] <- c("less","greater")[which(combinedpvaluelog2 == min(combinedpvaluelog2,na.rm=T))[1]]
          collationpvaluecountlog2[i,1] <- pvaluescountlog2[which(combinedpvaluelog2 == min(combinedpvaluelog2,na.rm=T))[1]]
          collationpvaluestartlog2[i,1] <- combinedpvaluestart[which(combinedpvaluelog2 == min(combinedpvaluelog2,na.rm=T))[1]]
          
          ####now combine consecutive significant ratios with same direction as indicated by selected combined pvalues 
          sum <- 0
          collationcounter <- 0
          counter <- 0
          direction <- 0 ### -1 if lower ratio < 0 and 1 if ratio > 0
          relevanttemps <- NA
          Values <- matrix(nrow = 0,ncol = 3) ### groups for data collation are summarized here to later find the optimum
          colnames(Values) <- c("log2diff","count","start")
          for(temp in cfg_info$temp_gradient) 
          {
            test <- subset(summary_data_combined[[ndata]], select = c(paste("meanlog2Ratio",temp,"C",sep=""),paste("pValue_less",temp,"C",sep=""),paste("pValue_greater",temp,"C",sep="")))[i,]
            colnames(test) <- c("log2diff","pValue_less","pValue_greater")
            
            if(!is.na(test[1]) && !is.na(test[2]) && test[2] <= cfg_info$pVal_collation_cutoff | !is.na(test[1]) && !is.na(test[3]) && test[3] <= cfg_info$pVal_collation_cutoff)
            {
              if(counter == 0) #Zähle wieviele Werte in einer Reihe zusammen unter dem Cutoff liegen
              {
                if(test[1] > 0){direction <- 1}
                if(test[1] < 0){direction <- -1}     
                collationcounter <- collationcounter + 1
                counter <- 1
                relevanttemps <- temp
                sum <- test[1]
              }else
              {
                if(direction == 1 && test[1]>0 | direction == -1 && test[1]<0)
                {
                  counter <- counter + 1
                  relevanttemps <- rbind(relevanttemps,which(cfg_info$temp_gradient == temp))
                  sum <- sum + test[1]
                }else ###significant difference but not in the correct direction then stop the combination here and start a new one from this point
                {
                  Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
                  colnames(Values) <- c("log2diff","count","start")
                  ######now go a step back and do a new collation at this point
                  if(test[1] > 0){direction <- 1}
                  if(test[1] < 0){direction <- -1}     
                  collationcounter <- collationcounter + 1
                  counter <- 1
                  relevanttemps <- temp
                  sum <- test[1]
                }
              }  
            }else if(counter != 0) ###if there was one value below the cutoff then save the pvalue and sumlog2diff of that collation and reset counter for next round of collation
            {
              Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
              colnames(Values) <- c("log2diff","count","start")
              counter <- 0
            }
            if(counter != 0 && which(cfg_info$temp_gradient == temp) == length(cfg_info$temp_gradient)) ###if also for the last temperature a significant difference was observed
            {
              Values <- rbind(Values,c(sum,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
              colnames(Values) <- c("log2diff","count","start")
              counter <- 0
            }
          }
          Values <- as.data.frame(Values)
          Values <- subset(Values,!is.na(log2diff))
          ##### select Group with abs(max log2ratio)
          selectedsum <- NA
          selectedcount <- NA
          selectedstart <- NA
          if(nrow(Values)>0)
          {
            max <- Values[which.max(Values$log2diff),]
            min <- Values[which.min(Values$log2diff),]
            if(abs(unlist(max$log2diff)) >= abs(unlist(min$log2diff)))
            {
              selectedsum <- unlist(max$log2diff)
              selectedcount <- unlist(max$count)
              selectedstart <- unlist(max$start)
            }else
            {
              selectedsum <- unlist(min$log2diff)
              selectedcount <- unlist(min$count)
              selectedstart <- unlist(min$start)
            } 
          }else ### if the pvalue cutoff was never reached, then use the value at most significant temp in the direction of selected collated pvalues
          {
            testdata <- summary_data_combined[[ndata]][i,which(grepl("meanlog2Ratio",colnames(summary_data_combined[[ndata]])))]
            testpvalues <- summary_data_combined[[ndata]][i,which(grepl(paste("pValue_",collationpvaluedirectionlog2[i,1],sep=""),colnames(summary_data_combined[[ndata]])))]
            selectedsum <- testdata[1,which(testpvalues == min(testpvalues,na.rm=T))]
            selectedcount <- 1
            selectedstart <- which(testpvalues == min(testpvalues,na.rm=T)) ###temp at max ratio
          }
          
          collationlog2ratio[i,1] <- selectedsum
          collationratiocount[i,1] <- selectedcount
          collationratiostart[i,1] <- selectedstart
          
          if(!is.na(selectedsum))
          {
            if(selectedsum > 0)
            {
              collationpValuelog2[i,1] <- combinedpvaluelog2[2]
            }else
            {
              collationpValuelog2[i,1] <- combinedpvaluelog2[1]
            }
          }
          
        }
        if(progressbar == T){setWinProgressBar(pb, i, label=paste( round(i/nrow(summary_data_combined[[ndata]])*100, 0)," % done (",i,"/",nrow(summary_data_combined[[ndata]]),")",sep = ""))}
      }
      if(progressbar == T){close(pb)}
      
      pdf(paste("Collated pvalue",ndata,".pdf",sep=""))
      plot(density(as.matrix(collationpValuelesslog2),na.rm=T),main="Collation pvalue_less")
      hist(as.numeric(as.matrix(collationpValuelesslog2)))
      plot(density(as.matrix(collationpValuegreaterlog2),na.rm=T),main="Collation pvalue_greater")
      hist(as.numeric(as.matrix(collationpValuegreaterlog2)))
      dev.off()
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationlog2ratio)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.log2Ratio"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpValuelog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpValuelesslog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue.less"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpValuegreaterlog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue.greater"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpvaluedirectionlog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue.direction"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpvaluecountlog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pvalue.count"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpvaluestartlog2)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue.start"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationratiocount)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.count"
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationratiostart)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.start"
      
      
    }
    
    if(cfg_info$pVal_collation_function == 0) ### collate significant
    {
      #######calculate significant summed log2 ratios and p-Values#########
      collationlog2ratio <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationpValue <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationcount <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      collationstart <- as.data.frame(matrix(ncol=1,nrow=nrow(summary_data_combined[[ndata]])))
      if(progressbar == T){pb <- winProgressBar(title = paste("Collating data:",ndata),label=paste( round(0/nrow(summary_data_combined[[ndata]])*100, 0),"% done"), min = 0,max = nrow(summary_data_combined[[ndata]]), width = 300)}
      for(i in 1:nrow(summary_data_combined[[ndata]])) 
      {
        sum <- 0
        pValue <- NA
        finalpValue <- NA
        collationcounter <- 0
        counter <- 0
        direction <- 0 ### -1 if lower ratio < 0 and 1 if ratio > 0
        relevanttemps <- NA
        Values <- matrix(nrow = 0,ncol = 4) ### groups for data collation are summarized here to later find the optimum
        colnames(Values) <- c("log2diff","pValue","count","start")
        for(temp in cfg_info$temp_gradient) 
        {
          test2 <- subset(summary_data_combined[[ndata]], select = c(paste("meanlog2Ratio",temp,"C",sep=""),paste("pValue_less",temp,"C",sep=""),paste("pValue_greater",temp,"C",sep="")))[i,]
          minindex <- 1+which(test2[1,2:3] == min(test2[1,2:3],na.rm=T))
          minindex <- ifelse(length(minindex)>0,minindex,2)
          test <- t(as.data.frame(c(test2[1,1],test2[1,minindex])))
          colnames(test) <- c("log2diff","pValue")
          rownames(test) <- c()
          
          if(!is.na(test[1]) && !is.na(test[2]) && test[2] <= cfg_info$pVal_collation_cutoff)
          {
            if(counter == 0) #Zähle wieviele Werte in einer Reihe zusammen unter dem Cutoff liegen
            {
              if(test[1] > 0){direction <- 1}
              if(test[1] < 0){direction <- -1}     
              collationcounter <- collationcounter + 1
              counter <- 1
              pValue <- test[2]
              relevanttemps <- temp
              sum <- test[1]
            }else
            {
              if(direction == 1 && test[1]>0 | direction == -1 && test[1]<0)
              {
                counter <- counter + 1
                pValue <- rbind(pValue,test[2])
                relevanttemps <- rbind(relevanttemps,temp)
                sum <- sum + test[1]
              }else ###significant difference but not in the correct direction then stop the combination here and start a new one from this point
              {
                if(counter > 1)
                {
                  #######combination of p-Values by brown´s method##############
                  #######collect respective rawdata#############################
                  datamatrix <- matrix(nrow = 0,ncol = 2*cfg_info$replicates)
                  for(relevanttemp in relevanttemps)
                  {
                    cond1cols <- which(grepl(paste(relevanttemp,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                    cond2cols <- which(grepl(paste(relevanttemp,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                    datamatrix <- rbind(datamatrix,as.numeric((summary_data_combined[[ndata]])[i,c(cond1cols,cond2cols)]))
                  }
                  #######calculate combined p-Value after brown´s method########
                  datamatrix <- as.data.frame(datamatrix[,colSums(is.na(datamatrix))==0])
                  if(ncol(datamatrix) > 1)
                  {
                    finalpValue <- empiricalBrownsMethod(data_matrix=datamatrix, p_values=pValue, extra_info=FALSE)
                  }else
                  {
                    finalpValue <- min(pValue)
                  }
                }else
                {
                  finalpValue <- pValue
                }
                Values <- rbind(Values,c(sum,finalpValue,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
                colnames(Values) <- c("log2diff","pValue","count","start")
                ######now go a step back and do a new collation at this point
                if(test[1] > 0){direction <- 1}
                if(test[1] < 0){direction <- -1}     
                collationcounter <- collationcounter + 1
                counter <- 1
                pValue <- test[2]
                relevanttemps <- temp
                sum <- test[1]
              }
            }  
          }else if(counter != 0) ###if there was one value below the cutoff then save the pvalue and sumlog2diff of that collation and reset counter for next round of collation
          {
            if(counter > 1)
            {
              #######combination of p-Values by brown´s method##############
              #######collect respective rawdata#############################
              datamatrix <- matrix(nrow = 0,ncol = 2*cfg_info$replicates)
              for(relevanttemp in relevanttemps)
              {
                cond1cols <- which(grepl(paste(relevanttemp,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                cond2cols <- which(grepl(paste(relevanttemp,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                datamatrix <- rbind(datamatrix,as.numeric((summary_data_combined[[ndata]])[i,c(cond1cols,cond2cols)]))
              }
              #######calculate combined p-Value after brown´s method########
              datamatrix <- as.data.frame(datamatrix[,colSums(is.na(datamatrix))==0])
              if(ncol(datamatrix) > 1)
              {
                finalpValue <- empiricalBrownsMethod(data_matrix=datamatrix, p_values=pValue, extra_info=FALSE)
              }else
              {
                finalpValue <- min(pValue)
              }
            }else
            {
              finalpValue <- pValue
            }
            
            Values <- rbind(Values,c(sum,finalpValue,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
            colnames(Values) <- c("log2diff","pValue","count","start")
            counter <- 0
          }
          if(counter != 0 && which(cfg_info$temp_gradient == temp) == length(cfg_info$temp_gradient)) ###if also for the last temperature a significant difference was observed
          {
            if(counter > 1)
            {
              #######combination of p-Values by brown´s method##############
              #######collect respective rawdata#############################
              datamatrix <- matrix(nrow = 0,ncol = 2*cfg_info$replicates)
              for(relevanttemp in relevanttemps)
              {
                cond1cols <- which(grepl(paste(relevanttemp,"C_cond1\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                cond2cols <- which(grepl(paste(relevanttemp,"C_cond2\\.",sep=""),colnames(summary_data_combined[[ndata]])))[1:cfg_info$replicates]
                datamatrix <- rbind(datamatrix,as.numeric((summary_data_combined[[ndata]])[i,c(cond1cols,cond2cols)]))
              }
              #######calculate combined p-Value after brown´s method########
              datamatrix <- as.data.frame(datamatrix[,colSums(is.na(datamatrix))==0])
              if(ncol(datamatrix) > 1)
              {
                finalpValue <- empiricalBrownsMethod(data_matrix=datamatrix, p_values=pValue, extra_info=FALSE)
              }else
              {
                finalpValue <- min(pValue)
              }
            }else
            {
              finalpValue <- pValue
            }
            
            Values <- rbind(Values,c(sum,finalpValue,counter,which(cfg_info$temp_gradient == relevanttemps[1])))
            colnames(Values) <- c("log2diff","pValue","count","start")
            counter <- 0
          }
        }
        Values <- as.data.frame(Values)
        
        Values$log2diff <- as.numeric(as.character(Values$log2diff))
        Values$pValue <- as.numeric(as.character(Values$pValue))
        Values$count <- as.numeric(as.character(Values$count))
        Values$start <- as.numeric(as.character(Values$start))
        
        Values <- subset(Values,!is.na(log2diff))
        
        ##### select Group with abs(max log2ratio)
        selectedsum <- NA
        selectedpVal <- NA
        selectedcount <- NA
        selectedstart <- NA
        if(nrow(Values)>0)
        {
          max <- Values[which.max(Values$log2diff),]
          min <- Values[which.min(Values$log2diff),]
          if(abs(unlist(max$log2diff)) >= abs(unlist(min$log2diff)))
          {
            selectedsum <- unlist(max$log2diff)
            selectedpVal <- unlist(max$pValue)
            selectedcount <- unlist(max$count)
            selectedstart <- unlist(max$start)
          }else
          {
            selectedsum <- unlist(min$log2diff)
            selectedpVal <- unlist(min$pValue)
            selectedcount <- unlist(min$count)
            selectedstart <- unlist(min$start)
          } 
        }else ### if the pvalue cutoff was never reached, then use the values at abs(maxlog2diff)
        {
          ####find temp with max meanlog2Ratio
          tmp <- summary_data_combined[[ndata]][i,which(grepl("meanlog2Ratio",colnames(summary_data_combined[[ndata]])))]
          tmptemp <- cfg_info$temp_gradient[which(abs(tmp) == max(abs(tmp),na.rm=T))] ###temp at max ratio
          if(length(tmptemp)>0)
          {
            selectedsum <- summary_data_combined[[ndata]][i,paste("meanlog2Ratio",tmptemp[1],"C",sep="")]
            selectedpVal <- min(c(summary_data_combined[[ndata]][i,paste("pValue_less",tmptemp[1],"C",sep="")],summary_data_combined[[ndata]][i,paste("pValue_greater",tmptemp[1],"C",sep="")]),na.rm=T)
            selectedcount <- 1
            selectedstart <- which(cfg_info$temp_gradient == tmptemp[1])
          }
          
        }
        if(!is.na(selectedpVal))
        {
          collationlog2ratio[i,1] <- selectedsum
          collationpValue[i,1] <- selectedpVal
          collationcount[i,1] <- selectedcount
          collationstart[i,1] <- selectedstart
        }else
        {
          collationlog2ratio[i,1] <- NA
          collationpValue[i,1] <- NA
          collationcount[i,1] <- NA
          collationstart[i,1] <- NA
        }
        if(progressbar == T){setWinProgressBar(pb, i, label=paste( round(i/nrow(summary_data_combined[[ndata]])*100, 0)," % done (",i,"/",nrow(summary_data_combined[[ndata]]),")",sep = ""))}
      }
      if(progressbar == T){close(pb)}
      
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationlog2ratio)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.log2Ratio"
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationpValue)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.pValue"
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationcount)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.count"
      summary_data_combined[[ndata]] <- cbind(summary_data_combined[[ndata]],collationstart)
      colnames(summary_data_combined[[ndata]])[length(colnames(summary_data_combined[[ndata]]))] <- "collation.ratio.start"
    }
  }
  return(summary_data_combined)
}