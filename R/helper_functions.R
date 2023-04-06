
allmaxDensity <- function(vals) ###returns x value of all maxima which have >= 80% value of total maxima
{
  density <- density(vals)
  ts_y<-ts(density$y)
  tp<-turnpoints(ts_y)
  maxima <- which(tp$peaks) ####all maxima
  totalmax <- max(density$y[maxima],na.rm=T)
  maxima <- maxima[which(density$y[maxima] >= 0.8*totalmax)]###only maxima with >=80% of total maxima
  return(density$x[maxima])
}

maxDensity=function(Data,from=NA,to=NA)
{
  if(all(Data==0,na.rm=T) | all(is.na(Data)))
  {
    return(0)
  }else
  {
    if(!is.na(from) && !is.na(to))
    {
      density <- density(Data,na.rm=T,from=from,to=to)
    }else
    {
      density <- density(Data,na.rm=T)
    }
    max <- density$x[which(density$y==max(density$y))]
    if(length(max)>1){max <- mean(max)}
    return(max)
  }
}

saveExcel = function(Data,Sheet_Names = NA,File)
{
  wb <- createWorkbook("Report")
  if(is.data.frame(Data) == 1)
  {
    if(!is.na(Sheet_Names))
    {
      addWorksheet(wb,Sheet_Names[1])
      writeData(wb, Sheet_Names[1],Data) 
    }else
    {
      addWorksheet(wb,"sheet_1")
      writeData(wb, "sheet_1",Data) 
    }
  }else
  {
    for(i in 1:length(Data))
    {
      if(any(is.na(Sheet_Names)))
      {
        addWorksheet(wb,paste("sheet_",i,sep=""))
        writeData(wb, paste("sheet_",i,sep=""),Data[[i]]) 
      }else
      {
        addWorksheet(wb,Sheet_Names[i])
        writeData(wb, Sheet_Names[i],Data[[i]]) 
      }
    }
  }
  
  saveWorkbook(wb, File, overwrite = TRUE)
}

removeWeakQuant = function(data_list_normalized,cfg_info) ###remove data with weak quantification
{
  POI <- subset(cfg_info$cfg_sheet$Proteins.of.interest,!is.na(cfg_info$cfg_sheet$Proteins.of.interest))
  for(i in 1:cfg_info$num_samples)
  {
    for(p in 1:nrow(data_list_normalized[[i]]))
    {
      #if(data_list_normalized[[i]]$gene_name[p] %not in% POI)
      {
        if(data_list_normalized[[i]]$qusm[p] < cfg_info$qusm_filter_plot | data_list_normalized[[i]]$qupm[p] < cfg_info$qupm_filter_plot)
        {
          data_list_normalized[[i]][p,which(grepl("rel_fc_protein",colnames(data_list_normalized[[i]])))] <- NA
          data_list_normalized[[i]][p,which(grepl("sumionarea_protein",colnames(data_list_normalized[[i]])))] <- NA
        }
      }
    }
  }
  return(data_list_normalized)
}

saveNormData = function(file,data_list,cfg_info,cond1_normalization_parameter_sheet,cond2_normalization_parameter_sheet) ###Save normalized data, each replicate and gradientpart separate
{
  wb <- createWorkbook("Normalized data")
  for(g in 1:cfg_info$num_gradientparts) #### go through data of specific gradient parts
  {
    dataindices <- which(cfg_info$cfg_sheet$Gradient == g)
    for(i in dataindices)
    {
      replicate <- which(dataindices == i)
      addWorksheet(wb, paste("Data.",replicate,".",g,sep=""))
      writeData(wb, paste("Data.",replicate,".",g,sep=""),data_list[[i]])
    }
  }
  addWorksheet(wb, "Norm_Cond1")
  addWorksheet(wb, "Norm_Cond2")
  writeData(wb, "Norm_Cond1",cond1_normalization_parameter_sheet)
  writeData(wb, "Norm_Cond2",cond2_normalization_parameter_sheet)
  saveWorkbook(wb, file, overwrite = TRUE)
}

replotAnalysis = function(path_to_analysis_cfg,path_to_output_folder) ###do ratio-based TPP data analysis
{
  ###################################################################################
  Analysisversion <- 1.0
  Internalizationcategories <<- c(0.25,0.5) #### Categories for internalization --> 0: < 0.1; 1: 0.1 < 0.25; 2: 0.25 < 0.5; 3: >= 0.5
  Thermalshiftcategories <<- c(1,2) #### Categories for thermalshift --> 0: < 0.5; 1: 0.5 < 1.0; 2: 1.0 < 2.0; 3: >= 2.0
  colors_melt_cond1 <<- c("red4","red3","red2","orangered4","orangered3","orangered2","indianred4","indianred3","indianred2","deeppink4","deeppink3","deeppink") ##in total colors for max 12 replicates
  colors_melt_cond2 <<- c("royalblue4","royalblue3","royalblue2","slateblue4","slateblue3","slateblue2","steelblue4","steelblue3","steelblue2","turquoise4","turquoise3","turquoise1") ##in total colors for max 12 replicates
  ###################################################################################
  
  ###Read AnalysisInfo-Sheet
  cfg_info <- readAnalysisSetup(path_to_analysis_cfg)
  
  setwd(path_to_output_folder)
  ###load results
  summary_data_combined <- list()
  summary_data_combined[[1]] <- read.xlsx(paste(path_to_output_folder,"\\Data summary - with abundance effect.xlsx",sep=""))
  summary_data_combined[[2]] <- read.xlsx(paste(path_to_output_folder,"\\Data summary - without abundance effect.xlsx",sep=""))
  norm_factors_meltingcurve <- list()
  if(cfg_info$normalize_melting_curve == 1)
  {
    temp <- read.xlsx(paste(path_to_output_folder,"\\Meltingcurve normalization - with abundance effect.xlsx",sep=""))
    temp <- temp[,-1]
    norm_factors_meltingcurve[["norm_factors_meltingcurve_ctrl_with_abundance"]] <- temp[1:cfg_info$replicates,]
    norm_factors_meltingcurve[["norm_factors_meltingcurve_treat_with_abundance"]] <- temp[(cfg_info$replicates+1):(2*cfg_info$replicates),]
    temp <- read.xlsx(paste(path_to_output_folder,"\\Meltingcurve normalization - without abundance effect.xlsx",sep=""))
    temp <- temp[,-1]
    norm_factors_meltingcurve[["norm_factors_meltingcurve_ctrl_without_abundance"]] <- temp[1:cfg_info$replicates,]
    norm_factors_meltingcurve[["norm_factors_meltingcurve_treat_without_abundance"]] <- temp[(cfg_info$replicates+1):(2*cfg_info$replicates),]
  }
  plotResults(summary_data_combined,norm_factors_meltingcurve,cfg_info,"both",Internalizationcategories,Thermalshiftcategories,colors_melt_cond1,colors_melt_cond2)
}