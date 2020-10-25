readAnalysisSetup = function(path_to_analysis_cfg) ###read analysis setup information
{
  configInfo <- list()
  configInfo[["directory"]] <- dirname(path_to_analysis_cfg)
  configInfo[["cfg_sheet"]] <- read.xlsx(path_to_analysis_cfg,1)
  
  ####determine number of temperatures per TMT experiment
  num_temps <- length(which(grepl("^Temp[0-9]$",colnames(configInfo$cfg_sheet))))
  
  ##### read-in experimental settings (Temperatures and Temp.order)
  configInfo[["temp_gradient"]] <- as.character(configInfo[["cfg_sheet"]]$Temperatures)
  configInfo[["temp_gradient"]] <- configInfo[["temp_gradient"]][!is.na(configInfo[["temp_gradient"]])]
  ###Temperatur order per TMT experiment
  configInfo[["num_gradientparts"]] <- length(unique(configInfo[["cfg_sheet"]]$Gradient[!is.na(configInfo[["cfg_sheet"]]$Gradient)]))
  configInfo[["temp_order"]] <- matrix(nrow=num_temps,ncol=configInfo[["num_gradientparts"]])
  for(i in 1:configInfo[["num_gradientparts"]])
  {
    temp <- subset(configInfo[["cfg_sheet"]],Gradient == i)[1,]
    configInfo[["temp_order"]][,i] <- t(temp[1,which(colnames(configInfo[["cfg_sheet"]]) %in% paste("Temp",1:num_temps,sep=""))])
  }
  configInfo[["temp_order"]] <- as.data.frame(configInfo[["temp_order"]])
  colnames(configInfo[["temp_order"]]) <- paste("Temp_Gradientpart_",unique(configInfo[["cfg_sheet"]]$Gradient[!is.na(configInfo[["cfg_sheet"]]$Gradient)]),sep="")
  configInfo[["num_samples"]] <- length(configInfo[["cfg_sheet"]]$MSExperimentIDs[!is.na(configInfo[["cfg_sheet"]]$MSExperimentIDs)])
  
  dataset_to_gradientpart <- as.data.frame(matrix(ncol=2,nrow=configInfo[["num_samples"]]))
  colnames(dataset_to_gradientpart) <- c("Gradientpart","Replicate")
  dataset_to_gradientpart$Gradientpart <- configInfo[["cfg_sheet"]]$Gradient[!is.na(configInfo[["cfg_sheet"]]$Gradient)]
  dataset_to_gradientpart$Replicate <- configInfo[["cfg_sheet"]]$Replicate[!is.na(configInfo[["cfg_sheet"]]$Replicate)]
  configInfo[["dataset_to_gradientpart"]] <- dataset_to_gradientpart
  
  ####TMTchannel per condition1/condition and temperature
  configInfo[["tmt_channel_per_exp_condition1"]] <- as.data.frame(matrix(ncol=configInfo[["num_samples"]],nrow=num_temps))
  configInfo[["tmt_channel_per_exp_condition2"]] <- configInfo[["tmt_channel_per_exp_condition1"]]
  for(i in 1:configInfo[["num_samples"]])
  {
    configInfo[["tmt_channel_per_exp_condition1"]][,i] <- as.data.frame(t(configInfo[["cfg_sheet"]][i,which(grepl("condition1_TMT_channel",colnames(configInfo[["cfg_sheet"]])))]))
    configInfo[["tmt_channel_per_exp_condition2"]][,i] <- as.data.frame(t(configInfo[["cfg_sheet"]][i,which(grepl("condition2_TMT_channel",colnames(configInfo[["cfg_sheet"]])))]))
  }
  colnames(configInfo[["tmt_channel_per_exp_condition1"]]) <- configInfo[["cfg_sheet"]]$MSExperimentIDs[!is.na(configInfo[["cfg_sheet"]]$MSExperimentIDs)]
  rownames(configInfo[["tmt_channel_per_exp_condition1"]]) <- paste("Temp",1:nrow(configInfo[["tmt_channel_per_exp_condition1"]]),sep="")
  colnames(configInfo[["tmt_channel_per_exp_condition2"]]) <- configInfo[["cfg_sheet"]]$MSExperimentIDs[!is.na(configInfo[["cfg_sheet"]]$MSExperimentIDs)]
  rownames(configInfo[["tmt_channel_per_exp_condition2"]]) <- paste("Temp",1:nrow(configInfo[["tmt_channel_per_exp_condition1"]]),sep="")
  
  configInfo[["replicates"]] <- configInfo$num_samples/configInfo$num_gradientparts
  configInfo[["subcellular_localization"]] <- configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "subcellular_localization")] #normalization and filtering of data based on subcellular localization
  configInfo[["norm_in_log2"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "norm_in_log2")]) #data normalization using log2 relfcs (1) or unlogged relfcs (0)
  configInfo[["normalize_melting_curve"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "normalize_melting_curve")]) #additional normalization of melting curves for melting curve fitting (1) or no additional normalization (0)
  configInfo[["qusm_filter_norm"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "qusm_filter_norm")]) #min qusm required to be used for norm.
  configInfo[["qupm_filter_norm"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "qupm_filter_norm")]) #min qupm required to be used for norm.
  configInfo[["relfc_cutoff_norm"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "relfc_cutoff_norm")]) #min relfc per temp which has to be exceed in cond1 or cond2 To be used for normalization
  configInfo[["remove_weak_quant"]] <- configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "remove_weak_quant")] #use also weakly quantified relcs for pvalue calculation (0) or remove these values (1)
  configInfo[["pVal_collation_cutoff"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "pVal_collation_cutoff")]) #pvalue cutoff for data collation
  configInfo[["plateau_ratio_cutoff"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "plateau_ratio_cutoff")]) #min ratio between plateaus of condition1 and condition2 to accept proteins with >= 50 % of abundance ratios significant in plateau range
  configInfo[["minstdev_factor_significance_function"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "minstdev_factor_significance_function")]) #multiplication factor for the stdev between replicates to define asymptotic convergence of the significance cutoff towards the y axis
  configInfo[["min_pvalue_significance_function"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "min_pvalue_significance_function")]) #minimal pvalue for the asymptotic convergence of the significance cutoff towards the x axis
  configInfo[["qusm_filter_plot"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "qusm_filter_plot")]) #min qusm required to be plotted
  configInfo[["qupm_filter_plot"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "qupm_filter_plot")]) #min qupm required to be plotted
  configInfo[["relfc_cutoff_plot"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "relfc_cutoff_plot")]) #min relfc per temp which has to be exceed in cond1 or cond2 To be plotted
  configInfo[["max_relfc_highesttemp_plot"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "max_relfc_highesttemp_plot")]) #max relfc at highest temp in one of both conditions to be plotted
  configInfo[["Rsq_cutoff"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "Rsq_cutoff")]) #min R² which has to be exceeded by melting curve fits for cond1 and cond2 at median 
  configInfo[["relnum_sigtemps_stdev_outlier"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "relnum_sigtemps_stdev_outlier")]) #maximal allowed relative number of significant temps with too high stdev compared to all proteins at respective temperature
  configInfo[["thermshift_cutoff"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "thermshift_cutoff")]) #min. thermal shift which has to be exceeded to be regarded as a significant thermal shift
  configInfo[["pvalue_thermshift_cutoff"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "pvalue_thermshift_cutoff")]) #pvalue cutoff for thermal shift to be regarded as a significant thermal shift
  configInfo[["abundance_cutoff"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "abundance_cutoff")]) #min relfc change which has to be exceeded to be regarded as a significant abundance change
  configInfo[["pvalue_abundance_cutoff"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "pvalue_abundance_cutoff")]) #pvalue cutoff for abundance changes to be regarded as a significant abundance change
  configInfo[["min_pvalue_BH_corrected"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "min_pvalue_BH_corrected")]) #min accepted pvalue after BH-correction of collated pvalues
  configInfo[["plateau_stdev_plot_filter"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "plateau_stdev_plot_filter")]) #plateau and Stdev filter for plot
  configInfo[["color_significant_by_shift_abundance"]] <- as.numeric(configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "color_significant_by_shift_abundance")]) #0 if significant proteins should be without any coloring 1 if proteins should be colored by thermal shift and abundance
  configInfo[["condition_name_1"]] <- configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "condition_name_1")] #Name for condition 1
  configInfo[["condition_name_2"]] <- configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "condition_name_2")] #Name for condition 2
  configInfo[["pVal_collation_function"]] <- configInfo$cfg_sheet$Filter_parameter_value[which(configInfo$cfg_sheet$Filter_parameter_name == "pVal_collation_function")] #how pvalues are collated; 0 = most significant; 1 = sliding window; 2 = all available
  
  
  return(configInfo)
}