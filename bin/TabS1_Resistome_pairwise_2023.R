#####
#Mikael Lenz Strube & Katrine Tams
#A script to make NMDS plots for each sampling date of the 16S data colored by treatment time and a PERMANOVA analysis of the same groups

#devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

rm(list = ls())
source("bin/functions.R")

library(vegan)
library(beeswarm)
library(RColorBrewer)
library(phyloseq)
library(pairwiseAdonis)


addEnvfit=F
writeTable=F

#Remove specific genes
remove_gene=TRUE
genes_to_remove=c("TolC1", "VatE", "water", "X16s", "X16s.1", "X16s.2") # TolC is a e.coli housekeeoing gene, not an ARG, and VatE had technical issues on one of the chips squewing week 10 

allData=CleanData_qPCR(PATH="data/fluidigm_runs/",removeNul = F,condenseB = T)

length(allData$FULL)
length(unique(allData$FULL))

#Alter sampling date for plotting purpose 
allData$DATE=gsub("2020-08-26","2020-08-25", allData$DATE)
allData$PIG_DATE=gsub("2020-08-26","2020-08-25",allData$PIG_DATE)
allData$FULL=as.vector(as.character(allData$FULL))

#Removing unvanted gene 
if(remove_gene){
  allData3=allData
  allData3$DATA=allData3$DATA[,-which(colnames(allData$DATA) %in% genes_to_remove)]
  allData=allData3
}

##Adding treatment data 
#Reading in meta data 
meta_data <- read.delim("data/meta_data.txt")
meta_data$DATE=as.Date(meta_data$DATE,format = "%d-%m-%Y" )
meta_data$pig_date=str_c(meta_data$PIG, '_', meta_data$DATE)


#Creating placeholder for metadata
allData$OUA=allData$FULL
allData$Stald=allData$FULL
allData$TreatTime=allData$FULL
allData$Week=allData$FULL
allData$BGROUP=allData$FULL
allData$Date_dif_0=allData$FULL
allData$SOW=allData$FULL
allData$Note=allData$FULL





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
#save(allData, file="data_clean/allData_Resistome")


#Ensuring the dates are in order
uniq_dates=as.character(unique(allData[["DATE"]]))
uniq_weeks=unique(allData$Week)[order(uniq_dates)]
uniq_dates=uniq_dates[order(uniq_dates)]

##Making vectors for collecting PERMANOVA data
Rsquared=c()
Rweeks=c()
Rdates=c()
pValue=c()
pAdjusted=c()
pdates=c()
pweeks=c()


#Load PERMANOVA values
#load("data_clean/PERMANOVA_RESISTOME_BGROUP_with_week_0")
#Or run the loop to create the PERMANOVA values

PW_DF=data.frame(DATE="",WEEK="",ADO="",
                 "CvsB_T"="","CvsB_U"="","CvsA"="",  "B_TvsB_U"="", "B_TvsA"="",  "B_UvsA"="",
                 "CvsB_TR"="","CvsB_UR"="","CvsAR"="",  "B_TvsB_UR"="", "B_TvsAR"="",  "B_UvsAR"="",
                 BDISP=""
)

date=uniq_dates[1]
for (date in uniq_dates){
  sample_date=allData$DATE[which(allData$DATE==date)[1]]
  sample_week=allData$Week[which(allData$DATE==date)[1]]
  
  print(sample_date)
  print(sample_week)
  
  #Index for samples at given date
  dateIndex=which(allData[["DATE"]]==date)
  
  #To ensure there is at least two groups to run the PERMANOVA
  #Creating dataframes for the given date
  dateData=allData[["DATA"]][dateIndex,]
  dateOUA=allData[["OUA"]][dateIndex]
  dateBGROUP=allData[["BGROUP"]][dateIndex]
  dateNote=allData[["Note"]][dateIndex]
  
  
  #Calculating PERMANOVA for given date
  ADO=adonis2(dateData~factor(dateBGROUP), parallel = 1, permutations = 2000) #
  
  PW1=pairwise.adonis(dateData,factor(dateBGROUP),perm = 1000)
  
  if(any(PW1$p.value<0.05 ) ) {
    PW1=pairwise.adonis(dateData,factor(dateBGROUP),perm = 2000)
  }
  
  
  
  #Calculating the betadispersion between the groups
  DISP=anova(betadisper(d = vegdist(dateData),group = factor(dateBGROUP)))
  #print(betadisper(d = vegdist(dateData),group = factor(dateBGROUP)))
  
  dummy=c(date,sample_week,ADO$`Pr(>F)`[1],
          PW1$p.adjusted,PW1$R2, DISP$`Pr(>F)`[1])
  
  PW_DF=rbind(PW_DF,dummy)
}

PW_DF_fixed=PW_DF

PW_DF_fixed[,-c(1,2)]=apply(PW_DF_fixed[,-c(1,2)],2, as.numeric)

PW_DF_fixed=PW_DF_fixed[-1,]

if(writeTable) {
  write.table(PW_DF_fixed,file = "Results/Table_S1.csv", sep=";", row.names = F)
}