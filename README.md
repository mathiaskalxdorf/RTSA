# RTSA - Ratio-based thermal shift assay analysis for cell surface thermal proteome profiling data

<p align="center"> 

<img src="images/CSTPP analysis workflow.png" style="width: 50%; height: 50%"/>

</p>

### Description

Ratio-based thermal shift assay (RTSA) enables a more robust and sensitive detection of effects on thermal stability and in parallel on protein abundance in thermal proteome profiling data. RTSA was developed for the cell surface thermal proteome profiling (CS-TPP, Kalxdorf et al., Nature Methods, 2020). 

### Installation

The following installations are required:

 

 - [R](https://cran.r-project.org/bin/windows/base/) (Version 4.02 or above)
 - [Rtools](https://cran.r-project.org/bin/windows/Rtools/history.html) (select required version)
 - Optional: [RStudio](https://rstudio.com/products/rstudio/download/)

During installation, please keep default settings.

Next, we install the RTSA package from GitHub


```r

install.packages("devtools")

devtools::install_github("mathiaskalxdorf/RTSA")

```

If everything wents fine, RTSA should be installed and is ready to be used.

### Usage example
An example data set is stored in the GitHub repository in the folder "Example data". Download and store this folder somewhere suitable.

With the following lines we load the RTSA package and start the analysis:


```r

library(RTSA)
runRTSAnalysis()
```

A file-selection window opens. Please navigate to the downloaded folder and select the "RTSA_Analysis_config.xlsx" file.
In this folder, the proteomics data should be stored as well. The file names should follow the following rules:
PFC_"MSExperimentID"_PF_"PF-Number"_QV_"QV-Number"_"PNAME".txt
Please check how this naming corresponds to columns in the "RTSA_Analysis_config.xlsx" file. Furthermore, please check which columns should be present in the data files for the analysis.

The analysis config file containes all configuration parameters required to run the analysis. It contains the following columns:

 - MSExperimentIDs: MS Experiment ID of corresponding proteomics data file
 - Comment: Optional. Can be used to add some general comment on the respective data
 - PName: Purification ID, can be a unique sample ID
 - PF: Protein flagger number (can be kept at 1)
 - QV: Quantification number (can be kept at 1)
 - TempX: Multiple columns from Temp1-X. Per TMT(pro)-Experiment we can cover up to 5 (8) individual temperatures for controls and treatments. Numerical indication per sample which temperature from column "Temperatures" is covered in respective temperature block. In this example, data in row 1 covers 5 temperatures, is the first (Gradient = 1) of two data files for replicate 1 (Replicate = 1) and covers the temperatures 1 = 37 °C, 2 = 44 °C, 3 = 47.6 °C, 8 = 62.2 °C, and 9 = 65.6 °C.
 - Gradient: To enable coverage of more than 5 (8) temperatures per replicate, in this column it can be indicated if current sample is one of several proteomics data for a replicate.
 - Replicate: Indicate which replicate the respective data represents
 - Temperatures: List of temperatures covered in total in this data set. Each row represents one covered temperature counted from 1 - X.
 - TempX_condition1_TMT_channel: Multiple columns for Temp1-X. Per row indicate which TMT channel represents the respective temperature block for control condition.
 - TempX_condition2_TMT_channel: Multiple columns for Temp1-X. Per row indicate which TMT channel represents the respective temperature block for treatment condition.
 - Proteins.of.interest: Optional. List proteins of interest which should not be excluded by any filtering and which should be highlighted in the final volcano plots. Each row should specify a protein of interest (gene name in the proteomics data tables).
 - Filter_parameter_name: List of several filter and config parameter names. Don´t change.
 - Filter_parameter_value: Corresponding values for respective filter and config parameters.
 - Filter_parameter_description: Description of respective filter and config parameters. Don´t change.

When analysis is done, a new subfolder will be created containing the analysis result files. The important result files are:
 
 - Data summary - with abundance effect.xlsx: Containes all analysis results (statistics, log2 fold-changes, passed filter criteria) per protein considering treatment-induced protein abundance effects.  
 - Data summary - without abundance effect.xlsx: Containes all analysis results (statistics, log2 fold-changes, passed filter criteria) per protein removing treatment-induced protein abundance effects.  
 - Induced abundance effect.pdf: Volcano plot showing proteins with treatment-induced abundance changes.
 - Result plots - with abundance effects.pdf: Volcano plot and reconstructed melting curves for significantly affected proteins considering treatment-induced abundance effects.
 - Result plots - without abundance effects.pdf: Volcano plot and reconstructed melting curves for significantly affected proteins removing treatment-induced abundance effects.
 - Significant proteins - with abundance effect.xlsx: Filtered data summary table for significantly affected proteins considering treatment-induced abundance effects.
 - Significant proteins - without abundance effect.xlsx: Filtered data summary table for significantly affected proteins removing treatment-induced abundance effects.
 - Significant regulated proteins.xlsx: Filtered data summary table for proteins which show significant treatment-induced abundance changes.