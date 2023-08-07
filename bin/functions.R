#Katrine Tams - modified from Mikael Lenz Strube and Henrik Spiegelhauer 
# 
#Turns SSI ASV table into a useable format
#
#Arguments
# PATH:       path of CSV data file
# PATH_meta   path of meta data xlsx file
# normalize:  normalize data into 100.000 reads per sample
# minMean:    cutoff for minimum mean abundance of ASVs
# ConRegEx:   regex matches to controls and other undesirables
#
#Returns a list with the data and the two vectors of meta-data
#This script returns only pig data, no sows


#Note for handling dates:
#   %d -> Day
#   %m -> numeric month
#   %b -> abbreviated months
#   %B -> full month name
#   %y -> 2-digit year
#   %Y -> 4-digit year


# rm(list = ls())
# knitr::opts_knit$set(root.dir = here::here())
# source("bin/functions.R")


# PATH="data_raw/Preliminary_tax_and_counts_seq44_45_46_47.csv"
# PATH_meta = "data_raw/meta_raw_2019.03.26 OUA data kohorte ORG.xlsx"
# normalize=F
# minMeans=-1
# ConRegEx="Neg|NK|NEC|Pos|Neg|7858.2"



# library(tidyverse)
# library(readxl)
# library(phyloseq)
# library(igraph)
# 
# # rm(list=ls())
# PATH = "data/Tax_and_counts_seq55_56.csv"
# normalize = T
# minMeans=-1
# ConRegEx="NEG|Neg|NC|NEC|Pos|Ossa"
# removeNul=F
# condenseB=F
# meta_path="data/meta_data.txt"
# use_meta=T #use meta to clean rather than from ASV table
# debug=T

CleanData16S=function(PATH, normalize=F,minMeans=-1, ConRegEx="Neg|NEG|NC|NEC|Pos|Ossa", removeNul=F, condenseB=F, meta_path="data/meta_data.txt", use_meta=T,debug=F ) {
  
  
  META=read.table(meta_path, sep="\t", header = T)
  
  #read data
  data0=read.table(PATH, sep=";", header = T)
  
  #finding taxa-cols
  whichIsTaxa=c(grep("Kingdom",colnames(data0)): NCOL(data0))
  
  #removing taxa and seq
  data1=data0[,-c(1,2,whichIsTaxa)]
  
  #saving taxa
  ASV0=data0[,whichIsTaxa]
  
  #saving seqs
  seqs=data0$ASV
  
  #making asvs for later
  ASV1=do.call(paste, c(ASV0, sep = "_"))
  
  #transposing data
  data2=data.frame(t(data1))
  
  #adding ASVs
  colnames(data2)=ASV1
  
  #Removing Negatives and ossa ConRegEx="Neg|NEG|NC|NEC|Pos|Ossa"
  negIndx=grep(ConRegEx,rownames(data2))
  data3=data2[-negIndx,]
  
  firstWrong=grep("25_08_2020",rownames(data3))
  rownames(data3)[firstWrong]=gsub("25_08_2020","25.08.2020_nul",rownames(data3)[firstWrong])
  
  secondWrong=grep("87_24_11_2020",rownames(data3))
  rownames(data3)[secondWrong]=gsub("87_24_11_2020","87_24.11.2020_C2",rownames(data3)[secondWrong])
  
  if(removeNul) {
    
    isNul=grep("nul", rownames(data3))
    data3=data3[-isNul,]
    
  }
  
  #find empty ASVs, e.g. ASVs only in controls
  emptyCols=which(colSums(data3)==0)
  data4=data3[,-emptyCols]
  
  if(minMeans>0) {
    lowCols=which(colMeans(data4)<minMeans)
    data4=data4[,-lowCols]
    cat(paste(length(lowCols)," ASVs removed.\n"))
  }
  
  #if normalize, otherwise return data4
  if(normalize){
    data_clean=t(apply(t(data4), 
                       FUN = function(x) round(100000 * x/sum(x), 0), MARGIN = 2))
  } else {
    data_clean=data4
  }
  
  
  #should also deal with meta in case
  if(condenseB) {
    rownames(data_clean)=gsub("B1.1","B1", rownames(data_clean))
    rownames(data_clean)=gsub("B1.2","B1", rownames(data_clean))
  }
  
  which(META$FULL=="X269_16.12.2020_B1.2_S155")
  which(rownames(data_clean)=="X269_16.12.2020_B1.2_S155")
  
  #double check that all meta exists in data
  META_check=data.frame(META=META[match( gsub("_v3v4","",rownames(data_clean)),META$FULL),4],
                        DATA=c(gsub("_v3v4","",rownames(data_clean))))
  
  #sort meta by data
  META_sorted=META[match( gsub("_v3v4","",rownames(data_clean)),META$FULL),]  
  
  META_sorted$ISNA=is.na(META_sorted$PIG)

  if(use_meta) {
    #work out what data is not in meta
    isNA_meta=which(is.na(META_sorted$OUA) | META_sorted$GROUP=="NA" )
    
    #remove it
    data_clean=data_clean[-isNA_meta,]
    META_sorted=META_sorted[-isNA_meta,]
    
    #Extract metadata from meta
    pig  =factor(gsub("X([0-9]*).([0-9\\.]*).([A-Za-z0-9\\.]*)_.*","\\1", META_sorted$PIG))
    date =factor(gsub("X([0-9]*).([0-9\\.]*).([A-Za-z0-9\\.]*)_.*","\\2", META_sorted$DATE))
    group=factor(gsub("X([0-9]*).([0-9\\.]*).([A-Za-z0-9\\.]*)_.*","\\3", META_sorted$GROUP_B))    
  } else {
    #Extract metadata from names
    pig  =factor(gsub("X([0-9]*)_([0-9\\.]*)_([A-Za-z0-9\\.]*)_.*","\\1", rownames(data_clean)))
    date =factor(gsub("X([0-9]*).([0-9\\.]*)_([A-Za-z0-9\\.]*)_.*","\\2", rownames(data_clean)))
    group_old=factor(gsub("X([0-9]*)_([0-9\\.]*)_([A-Za-z0-9\\.]*)_.*","\\3", rownames(data_clean)))
  }
  
  #date is ok
  if(debug) print(table(date, useNA = "ifany"))
  
  
  # table(group_old)
  # table(group)
  
  if(!use_meta) {
    #Corecting typo
    index_to_correct=which(pig=="235" & group=="C2")
    group[index_to_correct]="A1"
    
    
    #Grouping C1 and C2 together 
    group=gsub(group, pattern = "C2", replacement = "C1")
    group=gsub(group, pattern = "C1", replacement = "C")
    
    #Renaming A1 and B1 into A and B 
    group=gsub(group, pattern = "A1", replacement = "A")
    group=gsub(group, pattern = "B1", replacement = "B")
  }
  
  #Creating unique samplesname 
  date=as.Date(as.character(date),format = "%d-%m-%y" )
  pig_date=str_c(pig, '_', date)
  
  #date is still ok
  if(debug) print(table(date, useNA = "ifany"))
  
  data_clean=as.data.frame(data_clean)
  
  returnList=    list(PIG=pig,
                      DATE=date,
                      PIG_DATE=pig_date,
                      GROUP=group,
                      FULL=rownames(data_clean),
                      DATA=data_clean)
  
  
  if(debug) print(table(returnList$DATE, useNA = "ifany"))
  
  #return list with data and sample names
  return(
    returnList
  )
}




meta_data_clean=function(PATH_meta="~/GitHub/Transfer_study/data/meta_data.txt", DATA) {
  allData=DATA
  #Alter sampling date for plotting purpose 
  allData$DATE=gsub("2020-08-26","2020-08-25", allData$DATE)
  allData$PIG_DATE=gsub("2020-08-26","2020-08-25",allData$PIG_DATE)
  
  ##Adding treatment data 
  #Reading in meta data 
  meta_data <- read.delim(PATH_meta)
  meta_data$DATE=as.Date(meta_data$DATE,format = "%d-%m-%Y" )
  meta_data$pig_date=str_c(meta_data$PIG, '_', meta_data$DATE)
  
  #Creating placeholder for metadata
  allData$OUA=allData$PIG
  allData$Stald=allData$PIG
  allData$TreatTime=allData$PIG
  allData$Week=allData$PIG
  allData$BGROUP=allData$PIG
  allData$Date_dif_0=allData$PIG
  allData$SOW=allData$PIG
  allData$Note=allData$PIG
  
  
  #Looping over data to insert metadata
  for(i in 1:length(allData$PIG_DATE)){
    if(allData$PIG_DATE[i] %in% meta_data$pig_date){
      metaIndex=which(meta_data$pig_date==allData$PIG_DATE[i])
      allData$OUA[i]=meta_data$OUA[metaIndex]
      allData$Stald[i]=meta_data$Stald[metaIndex]
      allData$TreatTime[i]=meta_data$Behandlingsdato[metaIndex]
      allData$Week[i]=meta_data$Week[metaIndex]
      allData$BGROUP[i]=meta_data$GROUP_B[metaIndex]
      allData$Date_dif_0[i]=meta_data$DATE_dif_nul[metaIndex]
      allData$SOW[i]=meta_data$So[metaIndex]
      allData$Note[i]=meta_data$Note[metaIndex]
      
      
    }else{
      allData$OUA[i]=NA
      allData$Stald[i]=NA
      allData$TreatTime[i]=NA
      allData$Week[i]=NA
      allData$BGROUP[i]=NA
      allData$Date_dif_0=NA
      allData$SOW=NA
      allData$Note=NA
      
      
    }
  }
  
  return(allData)
  
}



# PATH_meta = "data_raw/meta_raw_2019.03.26 OUA data kohorte ORG.xlsx"
# pigs_for_grouping = COOC_data$PIG_DATE
# pigs_for_grouping = Data_qPCR$sample

GroupTreatment_meta=function(PATH_meta, pigs_for_grouping) {
  ##Cleaning meta data
  
  #Handling data in different excell sheets
  meta_data_sheets <- excel_sheets(PATH_meta)
  
  
  df_orange <- read_excel(path = PATH_meta, sheet="Grise_orange")
  df_behandlede <- read_excel(path = PATH_meta, sheet="Grise_behandlede")
  
  
  #Reading in relevant meta data: OUA status, pig id and date
  df_orange <- df_orange %>%
    select(id_gris,dato, oua_2) %>%
    rename(oua=oua_2)
  
  df_behandlede <- df_behandlede %>%
    select(id_gris,dato, oua)
  
  
  #Binding the two dataset from the two different excell sheets together
  meta_data <- bind_rows(df_orange, df_behandlede) %>%
    unite("sample",id_gris,dato, sep="-", remove=FALSE) %>%
    rename("id_pig"=id_gris,"date"=dato)
  
  #Removing duplikate (sample from df_behandlede)
  meta_data <- distinct(meta_data,sample, .keep_all = TRUE)
  
  #Keeping track of removed data
  removed<- subset(meta_data, meta_data$oua == 3 | meta_data$oua == 4)
  
  #Removing breeding pigs OUA=3
  meta_data<-subset(meta_data, meta_data$oua != 3)
  
  #Removing dead pigs OUA=4
  meta_data<-subset(meta_data, meta_data$oua != 4)
  
  #Streamlining dates
  meta_data$date =as.character(meta_data$date)
  
  
  #Uniq pigs for lopp
  unqPIGS=unique(meta_data$id_pig)
  
  treatment_group_frame<- as.data.frame(cbind(meta_data, GROUP=meta_data$oua, Newlytreated=meta_data$oua))
  names(treatment_group_frame)<- c("PIG_DATE","PIG" , "DATE" ,"OUA","GROUP","Newlytreated")
  treatment_group_frame$GROUP<- NA
  treatment_group_frame$Newlytreated<- NA
  
  pig=145
  
  for( pig in unqPIGS){
    pigIndex=which(meta_data$id_pig==pig)
    pigFrame=cbind(treatment_group_frame[pigIndex,],Index=pigIndex)
    pigFrame <-pigFrame[order(pigFrame$DATE),]
    UnTreatFlag=TRUE
    
    if (all(pigFrame$OUA==1)){
      for (i in pigFrame$Index){
        treatment_group_frame$GROUP[as.numeric(i)]="Untreated"
        treatment_group_frame$Newlytreated[as.numeric(i)]="Untreated"
      }
    }
    
    if (all(pigFrame$OUA==0)){
      for (i in pigFrame$Index){
        treatment_group_frame$GROUP[as.numeric(i)]=pigFrame$DATE[1]
        TREATTIME=pigFrame$DATE[1]
        if (TREATTIME==treatment_group_frame$DATE[i]){
          treatment_group_frame$Newlytreated[as.numeric(i)]="Newlytreated"          
        } else{
          treatment_group_frame$Newlytreated[as.numeric(i)]="Treated"  
        }
      }
    }
    
    if (all(pigFrame$OUA==0)==FALSE  && all(pigFrame$OUA==1)==FALSE){
      for (i in 1:length(pigIndex)){
        if(pigFrame$OUA[i]==1){
          treatment_group_frame$GROUP[as.numeric(pigFrame$Index[i])]="Untreated"
          treatment_group_frame$Newlytreated[as.numeric(pigFrame$Index[i])]="Untreated"
        } else if (pigFrame$OUA[i]==0 && UnTreatFlag==TRUE){
          TREATTIME=pigFrame$DATE[i]
          treatment_group_frame$GROUP[as.numeric(pigFrame$Index[i])]=TREATTIME
          treatment_group_frame$Newlytreated[as.numeric(pigFrame$Index[i])]="Newlytreated"
          UnTreatFlag=FALSE
        } else if (pigFrame$OUA[i]==0 && UnTreatFlag==FALSE ){
          treatment_group_frame$GROUP[as.numeric(pigFrame$Index[i])]=TREATTIME
          treatment_group_frame$Newlytreated[as.numeric(pigFrame$Index[i])]="Treated"
        }
        
      }
    }
    
  }
  
  Groups <- treatment_group_frame[treatment_group_frame$PIG_DATE %in% pigs_for_grouping, ]
  Groups=Groups %>% arrange(factor(PIG_DATE, levels = pigs_for_grouping))
  
  #Making group names in weeks old piglets
  Groups$WEEK=Groups$GROUP
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-04-20", "Treated at 1 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-04-25", "Treated at 2 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-05-03", "Treated at 3 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-05-07", "Treated at 4 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-05-15", "Treated at 5 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-05-24", "Treated at 6 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-05-29", "Treated at 7 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-06-07", "Treated at 8 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-06-20", "Treated at 12 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-07-05", "Treated at 14 weeks", Groups$WEEK)
  # Groups$WEEK=ifelse(Groups$WEEK == "2018-09-20", "Treated at 26 weeks", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-04-20", "T1", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-04-25", "T2", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-05-03", "T3", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-05-07", "T4", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-05-15", "T5", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-05-24", "T6", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-05-29", "T7", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-06-07", "T8", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-06-20", "T12", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-07-05", "T14", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "2018-09-20", "T26", Groups$WEEK)
  Groups$WEEK=ifelse(Groups$WEEK == "Untreated", "UT", Groups$WEEK)
  
  #Making sample datenames into weeks old piglets  
  Groups$SAMPLEWEEK=Groups$DATE
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-04-20", "Week 1", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-04-25", "Week 2", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-05-03", "Week 3", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-05-07", "Week 4", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-05-15", "Week 5", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-05-24", "Week 6", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-05-29", "Week 7", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-06-07", "Week 8", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-06-20", "Week 12", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-07-05", "Week 14", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-09-06", "Week 24", Groups$SAMPLEWEEK)
  Groups$SAMPLEWEEK=ifelse(Groups$SAMPLEWEEK == "2018-09-20", "Week 26", Groups$SAMPLEWEEK)
  
  
  
  return(
    Groups
    
  )
}



SpecifyGroups=function(Data_input, OnlyMainDates=FALSE, SmallGroupsRemoved=FALSE, OnlyFirstDate=FALSE){
  allData=Data_input
  if(OnlyMainDates){
    allDataDATA=allData$DATA
    IndexDATA=which(names(allData)=="DATA")
    allDataMeta=allData[-IndexDATA]
    
    MainDates=c("2018-04-25", "2018-05-07", "2018-05-15", "2018-05-24",  "2018-05-29",  "2018-06-07",  "2018-06-20",  "2018-07-05")
    
    
    DateIndexToKeep=which(allData$DATE %in% MainDates)
    
    allDataDATA2=allDataDATA[DateIndexToKeep,]
    allDataMeta2=lapply(allDataMeta, FUN=function(x) x[DateIndexToKeep])
    
    allData3=append(list(DATA=allDataDATA2),allDataMeta2)
    allData=allData3
    
  }
  
  
  if(SmallGroupsRemoved){
    allDataDATA=allData$DATA
    allDataMeta=allData[-c(which(names(allData)=="DATA"))]
    
    UnWhantedGroupes=c("2018-05-03", "2018-05-24","2018-06-07","2018-06-20","2018-07-05","2018-09-20")
    
    
    PigIndexToRemove=which(allData$GROUP %in% UnWhantedGroupes)
    
    allDataDATA2=allDataDATA[-PigIndexToRemove,]
    allDataMeta2=lapply(allDataMeta, FUN=function(x) x[-PigIndexToRemove])
    
    allData3=append(list(DATA=allDataDATA2),allDataMeta2)
    allData=allData3
    
  }
  if(OnlyFirstDate){
    allDataDATA=allData$DATA
    allDataMeta=allData[-c(which(names(allData)=="DATA"))]
    
    FirstDate=c("2018-04-25")
    
    PigIndexToKeep=which(allData$DATE %in% FirstDate)
    
    allDataDATA2=allDataDATA[PigIndexToKeep,]
    allDataMeta2=lapply(allDataMeta, FUN=function(x) x[PigIndexToKeep])
    
    allData3=append(list(DATA=allDataDATA2),allDataMeta2)
    allData=allData3
    
  }
  
  return(allData)
  
}



#####
#Condences phylogeny to genus and returns a phyloseq object
#####

MakePhyloGenus=function(dateData,meta ) {
  tax_mat=str_split_fixed(colnames(dateData), pattern = "_",n = 8)
  rownames(tax_mat) = colnames(dateData)
  colnames(tax_mat) = c("Kingdom", "Phylum", "Class", "Order",
                        "Family", "Genus", "Species","Strain")
  
  OTU = otu_table(t(dateData), taxa_are_rows = TRUE)
  
  rownames(meta)=rownames(dateData)
  
  META = sample_data(meta)
  TAX = tax_table(tax_mat)
  physeq = phyloseq(OTU, META, TAX)
  
  #Species=tax_glom(physeq = physeq,taxrank = "Species")
  Genus  =tax_glom(physeq = physeq,taxrank = "Genus")
  #Family =tax_glom(physeq = physeq,taxrank = "Family")
  
  # rownames(Species@otu_table)=
  #   apply(tax_table(Species)[,-c(7:8)],MARGIN = 1,function(x) paste(x,collapse="_"))
  
  rownames(Genus@otu_table)=
    apply(tax_table(Genus)[,-c(7:8)],MARGIN = 1,function(x) paste(x,collapse="_"))
  
  # rownames(Family@otu_table)=
  #   apply(tax_table(Family)[,-c(6:8)],MARGIN = 1,function(x) paste(x,collapse="_"))
  return(Genus)
}




#-----------------------------------------------------
#Groups replicates and calc the median
# group_df must be an excel with two cols [group]defines name of replicates
# members defines the genes included in the group, must be a string seperated only by "," (no spaces)

###############################
#   group members
#   gene1 gene1.1,gene1.2,...  
#   ...   ...
#   geneN ...

###############################

#'data/working_data/GeneDuplicates.xlsx'

group_replicates <- function(df, group_path, drop=TRUE){
  
  group_df <- read_excel(group_path) %>%
    mutate(members = str_split(.$members, ','))
  print.data.frame(group_df, right =F)
  
  groups = group_df$group
  g=groups[1]
  
  for (g in groups){
    
    members <- unlist(group_df$members[group_df$group==g])
    df_members <- df[members]
    #median
    #df[g] <- rowSums(select(df,any_of(members)) / length(members))
    df[g] <- apply(df_members, 1, median, na.rm=T)
    
    
    if (drop){
      df <- select(df, -all_of(members))
    }
    
  }
  return(df)
}


# Create Id col and date col and revert dates.
date_revert <- function(x){
  str_list = str_split(x, "-", simplify = T)
  #print(str_list)
  if (length(str_list)==3){
    #print(str_list)
    out = paste(str_list[1], substring(str_list[3],3,6), substring(str_list[3],1,2), str_list[2], sep="-")
  }
  else {
    out = paste(str_list[1], str_list[4], str_list[3], str_list[2], sep="-")
  }
  return(out)
}



#Cleaning qPCR data and coresponding metadata

# PATH="data/fluidigm_runs/"
# ConRegEx="kohort|faecalis|oolmix"
# removeNul=F
# condenseB=T


CleanData_qPCR <- function(PATH, ConRegEx="kohort|faecalis|oolmix", removeNul=F, condenseB=T){
  
  library(tidyverse)
  library(readxl)
  
  
  #Cleaning fluidigm data
  
  files <- list.files(path = PATH, pattern = ".csv")
  # #files <- files[order(files,decreasing = T)]
  # #read data    
  # #starts at 11, read onlu 96 forward and keep headers
  # qPCR_data=read.csv(paste(PATH,files[1], sep=""), skip = 11, nrows = 96, stringsAsFactors = F, header = T)[,-c(1,99,100)]
  # qPCR_data$PLATEnr=-1
  # qPCR_data=qPCR_data[-c(1:NROW(qPCR_data)),]
  # colnames(qPCR_data)[1]="sample"
  # 
  # file=files[which(files=="Run6_only_vatE.csv")]
  # 
  # data=read.csv(paste(PATH,file, sep=""), skip = 11, nrows = 96, stringsAsFactors = F, header = T)[,-c(1,99,100)]
  # 
  # #make data numeric
  # data1=data.frame(sample=data[,1],apply(data[,-1],2, as.numeric))
  # VatE=as.vector(data1$VatE.2)
  # backup_names=as.vector(data1$sample)
  
  files=files[-which(files=="Run6_only_vatE.csv")]
  
  qPCR_data=c()
  
  for(file in files){
    #print(file)
    
    # if(file=="Run6_no_vatE.csv"){
    #   data=read.csv(paste(PATH,file, sep=""), skip = 11, nrows = 96, stringsAsFactors = F, header = T)[,-c(1,99,100)]
    #   
    #   #make data numeric
    #   data1=data.frame(sample=data[,1],apply(data[,-1],2, as.numeric))
    #   data1$PLATEnr=file
    #   #backup_names==data1$sample
    #   
    #   #print(row(data1))
    #   data1$VatE=VatE
    #   qPCR_data=rbind(qPCR_data, data1)
    #   
    # }else 
    
    
    
    data=read.csv(paste(PATH,file, sep=""), skip = 11, nrows = 96, stringsAsFactors = F, header = T)[,-c(1,99,100)]
    
    #make data numeric
    data1=data.frame(sample=data[,1],apply(data[,-1],2, as.numeric))
    data1$PLATEnr=file
    #print(row(data1))
    qPCR_data=rbind(qPCR_data, data1)
  }
  
  
  #Removing redundant gene
  redundantgenes=c("ermB.jucl.Gun","ermB.jucl.Gun.1","ermB.jucl.Gun.2")
  redundantgenesIndex= which(colnames(qPCR_data) %in% redundantgenes)
  qPCR_data=qPCR_data[,-redundantgenesIndex]
  
  #Grouping replikate genes
  #Gene duplikates are manually typed into GeneDuplicates.xlsx
  qPCR_data <- group_replicates(df=qPCR_data, group_path='data/GeneReplikates.xlsx', drop=T)
  
  
  #Change Sample format so it is equal to metadata format.
  
  
  #Remove rows with 16S control negative 
  
  # negative_x16 <- qPCR_data %>%
  #   filter(X16s > 900 | X16s.1 > 900 | X16s.2 > 900) %>% # NOTE: Burde det være OR og ikke and???
  #   select(sample)
  # 
  # qPCR_data <- qPCR_data %>%
  #   anti_join(negative_x16, by = "sample")
  
  negative_x16<-subset(qPCR_data, qPCR_data$X16s > 900 | qPCR_data$X16s.1 > 900 | qPCR_data$X16s.2 > 900) 
  qPCR_data<-subset(qPCR_data, qPCR_data$X16s < 900 & qPCR_data$X16s.1 < 900 & qPCR_data$X16s.2 < 900) 
  
  # Calc median and sd for X16 genes and remove outliers
  
  qPCR_data <- qPCR_data %>%
    mutate("avg" = (.$X16s +.$X16s.1 + .$X16s.2 )/3) %>%
    mutate("sd" = apply(select(.,contains("X16")), 1, sd)) %>%
    mutate("cv" = 100*sd/avg)
  
  df16S=cbind(qPCR_data$X16s, qPCR_data$X16s.1, qPCR_data$X16s.2)
  
  qPCR_data$median = apply(df16S, 1, median, na.rm=T)
  
  # pal=brewer.pal(12,"Paired")
  # pal[13]="black"
  # pal[14]="green"
  # 
  # cols=pal[factor(qPCR_data$PLATEnr)]
  # 
  # plot(qPCR_data$cv, col=cols, pch = 16)
  # legend(x = "topright",legend = levels(factor(qPCR_data$PLATEnr)), col = pal,pch = 16)
  # abline(h=5)
  # 
  # plot(qPCR_data$sd, col=cols, pch = 16)
  # legend(x = "topright",legend = levels(factor(qPCR_data$PLATEnr)), col = pal,pch = 16)
  # plot(qPCR_data$avg, col=cols, pch = 16)
  # legend(x = "topright",legend = levels(factor(qPCR_data$PLATEnr)), col = pal,pch = 16)
  
  
  
  outliers_16s <- qPCR_data %>%
    filter(cv > 5)
  
  
  
  qPCR_data <- anti_join(qPCR_data, outliers_16s, by="sample")
  
  # Normalise qPCR by median X16 genes
  qPCR_data <- qPCR_data %>%
    mutate(across(-c(sample, PLATEnr, median, sd, cv), function(x){2^(median-x)}))  %>% #Relative abundance of genes related to 16s 
    mutate(across(-c(sample, PLATEnr, median, sd, cv), function(x){ifelse(x<.0001,0,x)})) #Set very small numbers to 0
  
  
  # Name cleanup.
  #Removing Negatives and ossa 
  #unique(qPCR_data$sample)
  negIndx=grep(ConRegEx,qPCR_data$sample)
  qPCR_data1=qPCR_data[-negIndx,]
  
  firstWrong=grep("25_08_2020",qPCR_data1$sample)
  qPCR_data1$sample[firstWrong]=gsub("25_08_2020","25.08.2020_nul",qPCR_data1$sample[firstWrong])
  
  secondWrong=grep("87_24_11_2020",qPCR_data1$sample)
  qPCR_data1$sample[secondWrong]=gsub("87_24_11_2020","87_24.11.2020_C2",qPCR_data1$sample[secondWrong])
  
  if(removeNul) {
    
    isNul=grep("nul", qPCR_data1$sample)
    qPCR_data1=qPCR_data1[-isNul,]
    
  }
  
  
  
  if(condenseB) {
    qPCR_data1$sample=gsub("B1.1","B1", qPCR_data1$sample)
    qPCR_data1$sample=gsub("B1.2","B1", qPCR_data1$sample)  
  }
  
  qPCR_data1$sample=gsub(" ","", qPCR_data1$sample)
  qPCR_data1$sample=gsub("\\.","-", qPCR_data1$sample)
  
  
  meta=as.data.frame(str_split_fixed(qPCR_data1$sample, "_", 3))
  colnames(meta) <- c('pig','date','group')
  
  index_to_correct=which(meta$pig=="235" & meta$group=="C2")
  meta$group[index_to_correct]="A1"
  
  
  #Grouping C1 and C2 together 
  meta$group=gsub(meta$group, pattern = "C2", replacement = "C1")
  meta$group=gsub(meta$group, pattern = "C1", replacement = "C")
  
  #Renaming A1 and B1 into A and B 
  meta$group=gsub(meta$group, pattern = "A1", replacement = "A")
  meta$group=gsub(meta$group, pattern = "B1", replacement = "B")
  full_names=qPCR_data1[,1]
  cols_to_remove=which(colnames(qPCR_data1) %in% c( "sample","avg","sd","cv","median", "PLATEnr"))
  
  qPCR_data1=qPCR_data1[,-cols_to_remove]
  
  #Creating unique samplesname 
  meta$date=as.Date(meta$date,format = "%d-%m-%y" )
  pig_date=str_c(meta$pig, '_', meta$date)
  
  result=list(PIG=meta$pig,
              PIG_DATE=pig_date,
              DATE=meta$date,
              GROUP=meta$group,
              FULL=full_names,
              DATA=qPCR_data1)
  return(result)
}

#Filtering qPCR and 16S data for downstream coocurance analysis. Removing samples which are not in both data sets and removing variables with a procentage of zeroes higher than the zero_pct_cutuf value
pre_COOC_filtering_and_climb_genus <- function(Data_16S, Data_qPCR, zero_pct_cutuf=0.98){
  
  #Streamlinening dataframes by mergeing 
  df_16S=as.data.frame(cbind(Data_16S[["PIG_DATE"]],Data_16S[["OUA"]],Data_16S[["DATE"]],Data_16S[["DATA"]]))
  names(df_16S)[1] <- "sample_ID"
  names(df_16S)[2] <- "OUA"
  names(df_16S)[3] <- "DATE"
  df_qPCR=as.data.frame(cbind(Data_qPCR[["PIG_DATE"]],Data_qPCR[["OUA"]],Data_qPCR[["DATE"]],Data_qPCR[["DATA"]]))
  names(df_qPCR)[1] <- "sample_ID"
  names(df_qPCR)[2] <- "OUA"
  names(df_qPCR)[3] <- "DATE"
  
  merge16s_qpcr=merge(x=df_16S, y = df_qPCR, by.x = "sample_ID", by.y = "sample_ID")
  
  #merge16s_qpcr=merge(x=df_16S, y = data.frame(df_qPCR,suge=df_qPCR$sample) , by.x = "sample_ID", by.y = "sample", )
  
  #Seperating 16S and qPCR data after merging 
  COOC_data_16S=merge16s_qpcr[,which(grepl("Bacteria|Archaea",colnames(merge16s_qpcr)))]  
  COOC_data_qPCR=merge16s_qpcr[,-which(grepl("Bacteria|Archaea|DATE|OUA|sample|water|16",colnames(merge16s_qpcr)))]
  sample_ID=merge16s_qpcr$sample_ID
  OUA=merge16s_qpcr$OUA.y
  DATE=merge16s_qpcr$DATE.y
  
  #Removing rows with more zeroes than given cutoff value after merging
  #qPCR
  zeroPerc=c()
  for(i in 1:NCOL(COOC_data_qPCR)) {
    zeroPerc=c(zeroPerc,(length(which(COOC_data_qPCR[,i]==0))/length(COOC_data_qPCR[,i])))
  }
  tooManyZeroesIndx=which(zeroPerc>zero_pct_cutuf)
  length(tooManyZeroesIndx)
  COOC_data_qPCR_empty=COOC_data_qPCR[,tooManyZeroesIndx]
  COOC_data_qPCR_empty_names <- names(COOC_data_qPCR_empty)
  COOC_data_qPCR_empty_names_and_pnct <- cbind(names(COOC_data_qPCR_empty), zeroPerc[tooManyZeroesIndx])
  COOC_data_qPCR=COOC_data_qPCR[,-c(tooManyZeroesIndx)]
  
  #16S
  
  #Climbing up to genus
  COOC_data_16S_phyloseq=MakePhyloGenus(COOC_data_16S[,-1], data.frame(COOC_data_16S[,1]))
  COOC_data_16S_genus=data.frame(sample_ID,t(otu_table(COOC_data_16S_phyloseq)))
  
  
  zeroPerc=c()
  for(i in 1:NCOL(COOC_data_16S_genus)) {
    zeroPerc=c(zeroPerc,(length(which(COOC_data_16S_genus[,i]==0))/length(COOC_data_16S_genus[,i])))
  }
  tooManyZeroesIndx=which(zeroPerc>zero_pct_cutuf)
  length(tooManyZeroesIndx)
  
  COOC_data_16S_empty=COOC_data_16S_genus[,tooManyZeroesIndx]
  COOC_data_16S_empty_names <- names(COOC_data_16S_empty)
  
  COOC_data_16S_genus=COOC_data_16S_genus[,-c(tooManyZeroesIndx)]
  
  
  
  return(
    list(PIG_DATE=sample_ID,
         OUA=OUA,
         DATE=DATE, 
         COOC_data_qPCR=COOC_data_qPCR,
         COOC_data_16S_genus=COOC_data_16S_genus[,-1]
    )
  ) 
  
  
}




color.gradient <- function(x, colors=c("grey90","red"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

prettyGraph=function(network) {
  
  l2=layout_nicely(network) 
  #l2=layout_with_lgl(network) 
  
  argNodes=which(nchar(names(V(network)))<20)
  bacNodes=which(nchar(names(V(network)))>20)
  
  plot(x = network,  
       axes=F,rescale=T,layout=l2, 
       # === vertex
       vertex.color = c(rep(1,length(argNodes)),rep(2,length(bacNodes))) ,   # Node color
       #vertex.frame.color = ifelse(translationTable2$genus=="Virion",NA, NA),# Node border color
       #vertex.shape=ifelse(translationTable2$genus=="Virion","circle", "circle"),                        # One of “none”, “circle”, “square”, “csquare”, “rectangle” “crectangle”, “vrectangle”, “pie”, “raster”, or “sphere”
       vertex.size=5,                                # Size of the node (default is 15)
       vertex.size2=NA,                              # The second size of the node (e.g. for a rectangle)
       
       # === vertex label
       vertex.label=c(names(V(network))[argNodes],str_split_fixed(names(V(network))[bacNodes],"_",6)[,6]),#translationTable2$nodeNames,                 # Character vector used to label the nodes
       vertex.label.color="black",
       vertex.label.font=1,                          # Font: 1 plain, 2 bold, 3, italic, 4 bold italic, 5 symbol
       vertex.label.cex=.8,                          # Font size (multiplication factor, device-dependent)
       vertex.label.dist=0,                          # Distance between the label and the vertex
       vertex.label.degree=pi/2 ,                    # The position of the label in relation to the vertex (use pi)
       
       edge.color=color.gradient(log(E(network)$weight+1),colors = c("grey90","yellow","red")),#"grey50",                           # Edge color
       edge.width=E(network)$weight,#edge.betweenness(network)*0.01,                                 # Edge width, defaults to 1
       edge.arrow.size=1,                            # Arrow size, defaults to 1
       edge.arrow.width=1,                           # Arrow width, defaults to 1
       edge.lty="solid",                             # Line type, could be 0 or “blank”, 1 or “solid”, 2 or “dashed”, 3 or “dotted”, 4 or “dotdash”, 5 or “longdash”, 6 or “twodash”
       edge.curved=0.3                               # Edge curvature, range 0-1 (FALSE sets it to 0, TRUE to 0.5)
       #axes=T
       #xlim=c(min(l),max(l)),ylim=c(min(l),max(l))#, asp = 0
  )
}


shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  reposition_legend(p, 'center', panel=names)
}


lastTaxa=function(GENUS) {
  sigGenus2=GENUS
  #Genus names for plotting
  taxa=str_split_fixed(sigGenus2,pattern = "_",6)
  for(i in 1:NROW(taxa)) {
    if(any(grepl("NA",taxa[i,]))){
      phylLevelnotNA=min(which(taxa[i,]=="NA"))-1
      sigGenus2[i]=taxa[i,phylLevelnotNA]
    } else {
      sigGenus2[i]=taxa[i,6]
    }
  }
  return(sigGenus2)
}


