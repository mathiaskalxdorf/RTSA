readData = function(cfg_info,uniprot_annotations,progressbar=T,remove_contaminants=T) ###read data from lims based on information from analysis setup sheet
{
  dataList <- list()
  max <- cfg_info$num_samples
  if(progressbar == T){pb <- tcltk::tkProgressBar(title = "Collecting data ",label=paste( round(0/max*100, 0),"% done"), min = 0,max = max, width = 300)}
  for(i in 1:cfg_info$num_samples)
  {
    printf("Read MSExperiment Data %d",cfg_info$cfg_sheet[i,1])
    if("PF" %in% colnames(cfg_info$cfg_sheet)){PF <- cfg_info$cfg_sheet$PF[i]}else{PF <- 1}
    if("QV" %in% colnames(cfg_info$cfg_sheet)){QV <- cfg_info$cfg_sheet$QV[i]}else{QV <- 1}
    path <- paste("PFC_",cfg_info$cfg_sheet[i,1],"_PF",PF,"_QV",QV,"_",cfg_info$cfg_sheet[i,3],".txt",sep="")
    Data<-read.csv(path,header = TRUE,sep="\t")
    if(remove_contaminants==T & any(colnames(Data) == "focus_group"))
    {
      Data <- Data[which(Data$focus_group %not in% c("contaminant","potential contaminant")),]
    }
    indx <- which(grepl("proteinscore|totalpsm|^qupm$|^qusm$|means2i|sumionarea_protein|^rel_fc_protein|^log2_rel_fc_protein|^ms1intensity$|^ms1seqs$|^ms1maxminusmin$",colnames(Data)))
    Data[indx] <- lapply(Data[indx], function(x) as.numeric(as.character(x)))
    colnames(Data)[1] <- "gene_name"
    colnames(Data)[indx] <- gsub("\\.","",colnames(Data)[indx])
    colnames(Data)[indx] <- gsub("N$","L",colnames(Data)[indx])
    colnames(Data)[indx] <- gsub("C$","H",colnames(Data)[indx])
    Data <- left_join(Data,uniprot_annotations, by=c("gene_name" = "Gene.Name"))
    if(cfg_info$subcellular_localization == "all"){dataList[[i]] <- Data}
    if(cfg_info$subcellular_localization == "surface"){dataList[[i]] <- subset(Data,Localization %in% c("Plasma membrane","Secreted"))}
    if(progressbar == T){tcltk::setTkProgressBar(pb, i, label=paste( round(i/max*100, 0)," % done (",i,"/",max,")",sep = ""))}
  }
  if(progressbar == T){close(pb)}
  return(dataList)
}