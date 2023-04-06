library(xlsx)
library(openxlsx)
library(dplyr)
library(ggplot2)
library(plyr)
library(reshape2)
library(randomcoloR)
library(matrixStats)
library(graphics)
library(pastecs)
library(calibrate)
library(Hmisc)
library(stringr)
library(shape)
library(RColorBrewer)
library(EmpiricalBrownsMethod)
library(PECA)
library(data.table)


printf = function(s, ...) cat(paste0(sprintf(s, ...)), '\n')

`%not in%` <- function (x, table) is.na(match(x, table, nomatch=NA_integer_))

runRTSAnalysis = function(path_to_analysis_cfg=NA,progressbar=T,remove_contaminants=T) ###do ratio-based TPP data analysis
{
  ###################################################################################
  Analysisversion <<- 1.0
  Internalizationcategories <<- c(0.25,0.5) #### Categories for internalization --> 0: < 0.1; 1: 0.1 < 0.25; 2: 0.25 < 0.5; 3: >= 0.5
  Thermalshiftcategories <<- c(1,2) #### Categories for thermalshift --> 0: < 0.5; 1: 0.5 < 1.0; 2: 1.0 < 2.0; 3: >= 2.0
  colors_melt_cond1 <<- c("red4","red3","red2","orangered4","orangered3","orangered2","indianred4","indianred3","indianred2","deeppink4","deeppink3","deeppink") ##in total colors for max 12 replicates
  colors_melt_cond2 <<- c("royalblue4","royalblue3","royalblue2","slateblue4","slateblue3","slateblue2","steelblue4","steelblue3","steelblue2","turquoise4","turquoise3","turquoise1") ##in total colors for max 12 replicates
  ###################################################################################
  
  ###check if path to the analysis config xls file was defined
  if(is.na(path_to_analysis_cfg))
  {
    path_to_analysis_cfg <- file.choose(new = T)
  }
  
  ###get identifier name for analysis run
  time <- as.character(Sys.time())
  time <- gsub(" |:|-","",time)
  analysis_file_name <- substr(path_to_analysis_cfg,unlist(gregexpr("\\\\|/",path_to_analysis_cfg))[length(unlist(gregexpr("\\\\|/",path_to_analysis_cfg)))]+1,regexpr(".xlsx",path_to_analysis_cfg)[1]-1)
  analysis_run_name <- paste(analysis_file_name,"_",time,sep="")
  
  ###Read AnalysisInfo-Sheet
  cfg_info <- readAnalysisSetup(path_to_analysis_cfg)
  
  ###Create outputfolder
  setwd(cfg_info$directory)
  dir.create(analysis_run_name, showWarnings = FALSE)
  
  ###Read uniprot annotations
  uniprot_annotations <- read.xlsx(file.path(system.file("extdata", package = "RTSA"), "Human Proteome with Localization and Membrane Interaction.xlsx"))
  
  ####Read data from lims
  data_list <- readData(cfg_info,uniprot_annotations,progressbar,remove_contaminants)
  
  ###Change working directory and create outputfolder
  setwd(paste(cfg_info$directory,"/",analysis_run_name,sep=""))
  
  ####Check if reference temperature (37°C) of Control condition is used as reference channel
  for(i in 1:length(data_list))
  {
    ref_channel <- cfg_info$tmt_channel_per_exp_condition1[1,i]
    if(any(na.omit(data_list[[i]][,paste("rel_fc_protein_",ref_channel,sep="")]) != 1))##relative fold changes calculated to wrong reference channel
    {
      temp_sumionareas <- data_list[[i]][,which(grepl("sumionarea",colnames(data_list[[i]])))]
      temp_sumionareas[temp_sumionareas==0] <- NA
      
      rel_fcs <- temp_sumionareas/temp_sumionareas[,paste("sumionarea_protein_",ref_channel,sep="")]
      data_list[[i]][,which(grepl("rel_fc_protein_",colnames(data_list[[i]])))] <- rel_fcs
    }
  }
  
  ###Normalization of all Condition1 to each other at each individual temperature based on rel.fcs
  output_normalization_cond1 <- dataNormalizationofcond1(cfg_info,data_list,cfg_info$norm_in_log2)
  
  ###Normalization of all Condition2 to condition1 at each individual temperature based on rel.fcs
  output_normalization_cond2 <- dataNormalizationofcond2(cfg_info,output_normalization_cond1$data_list_norm,cfg_info$norm_in_log2)
  data_list_normalized <- output_normalization_cond2$data_list_norm
  
  ###Remove data with weak quantification
  if(cfg_info$remove_weak_quant == 1){data_list_normalized <- removeWeakQuant(data_list_normalized,cfg_info)}
  
  ###save data - normalized data with abundance effect
  saveNormData(file = "Normalized data - with abundance effect.xlsx",data_list_normalized,cfg_info,output_normalization_cond1$cond1_normalization_parameter_sheet,output_normalization_cond2$cond2_normalization_parameter_sheet)
  
  ###save data with/without abundance effect
  summary_data_list_normalized <- list() 
  summary_data_list_normalized[["with_abundance_effect"]] <- data_list_normalized
  
  ###calculate relfcs of cond2 relative to cond2 at reference temperature
  num_temps <- length(which(grepl("^Temp[0-9]$",colnames(cfg_info$cfg_sheet))))
  for(i in 1:cfg_info$num_samples)
  {
    colindexref <- which(colnames(data_list_normalized[[i]]) == paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[1,i],"_norm",sep=""))
    for(j in 2:num_temps)
    {
      colindextemp <- which(colnames(data_list_normalized[[i]]) == paste("rel_fc_protein_",cfg_info$tmt_channel_per_exp_condition2[j,i],"_norm",sep=""))
      data_list_normalized[[i]][,colindextemp] <- data_list_normalized[[i]][,colindextemp]/data_list_normalized[[i]][,colindexref]
    }
    data_list_normalized[[i]][,colindexref] <- data_list_normalized[[i]][,colindexref]/data_list_normalized[[i]][,colindexref]
  }
  
  ###save data - normalized data without abundance effect
  saveNormData(file = "Normalized data - without abundance effect.xlsx",data_list_normalized,cfg_info,output_normalization_cond1$cond1_normalization_parameter_sheet,output_normalization_cond2$cond2_normalization_parameter_sheet)
  summary_data_list_normalized[["without_abundance_effect"]] <- data_list_normalized
  
  ###Calculate ratios between condition2 and condition1, define if calculation is done for data with/without abundance effect or both
  summary_data_list_normalized <- calculateRatios(summary_data_list_normalized,cfg_info,fordata = "both") ###define if calculation is done for data with/without abundance effect or both using fordata which can be 1: without 2: with or 3: both
  
  ###combine datasets to single dataframe and calculate mean ratios per temperature
  summary_data_combined <- combineData(summary_data_list_normalized,cfg_info,"both")
  
  ###check if one temperature is consistantly missing e.g. due to sample loss and hence reporter ion channels are always empty
  exclude_temp <- NULL
  for(t in cfg_info$temp_gradient)
  {
    sel_cols <- which(grepl(paste(t,"C_cond",sep=""),colnames(summary_data_combined$with_abundance_effect)))
    
    temp_dat <- summary_data_combined$with_abundance_effect[,sel_cols]
    if(!any(!is.na(temp_dat)))##all missing
    {
      exclude_temp <- append(exclude_temp,t)
    }
  }
  if(!is.null(exclude_temp)) ###
  {
    #remove respective cols from combined tables
    for(t in exclude_temp)
    {
      sel_cols <- which(grepl(t,colnames(summary_data_combined$with_abundance_effect)))
      summary_data_combined[[1]] <- summary_data_combined[[1]][-sel_cols]
      summary_data_combined[[2]] <- summary_data_combined[[2]][-sel_cols]
    }
    cfg_info$temp_gradient <- cfg_info$temp_gradient[-which(cfg_info$temp_gradient %in% exclude_temp)]
  }
  
  ###calculate pValues for log2 relfcs at each temperature
  summary_data_combined <- calculatepvalues(summary_data_combined,cfg_info,"both",progressbar)
  
  ###collation of data points for temperatures
  summary_data_combined <- dataCollation(summary_data_combined,cfg_info,"both",progressbar)
  
  ###fit melting curves, calculate melting points and rsqares as well as determine melting plateaus
  output_meltingcurve_fitting <- meltingCurveFitting(summary_data_combined,cfg_info,"both",progressbar)
  summary_data_combined <- output_meltingcurve_fitting$summary_data_combined
  
  ###prepare filter criteria columns for data filtering
  summary_data_combined <- prepareDataFiltering(summary_data_combined,cfg_info,"both")
  
  ###save results
  Data <- list()
  Data[[1]] <- summary_data_combined[[1]]
  Data[[2]] <- cfg_info$cfg_sheet
  saveExcel(Data = Data,Sheet_Names = c("Data","Cfg_File"),File = "Data summary - with abundance effect.xlsx")
  Data[[1]] <- summary_data_combined[[2]]
  Data[[2]] <- cfg_info$cfg_sheet
  saveExcel(Data = Data,Sheet_Names = c("Data","Cfg_File"),File = "Data summary - without abundance effect.xlsx")
  
  ###plot results
  plotResults(summary_data_combined,output_meltingcurve_fitting$norm_factors_meltingcurve,cfg_info,"both",Internalizationcategories,Thermalshiftcategories,colors_melt_cond1,colors_melt_cond2)
  
}



