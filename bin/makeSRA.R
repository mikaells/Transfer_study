
#####
#Mikael Lenz Strube & Katrine Tams
#A script to make NMDS plots for each sampling date of the 16S data colored by treatment time and a PERMANOVA analysis of the same groups

rm(list = ls())
source("bin/functions.R")

library(vegan)
library(beeswarm)
library(RColorBrewer)
library(phyloseq)
library(stringr)

addEnvfit=F
makePNG=T
#namingScheme=c("A","B_T","B_U","C")
namingScheme=c("Treated","Mix-Treated","Mix-Untreated","Untreated")


allData=CleanData16S(PATH = "data/Tax_and_counts_seq55_56.csv",normalize = T,ConRegEx="NEG|Neg|NC|NEC|Pos",removeNul = F,condenseB = T)

unique(allData$PIG)

#Alter sampling date for plotting purpose 
allData$DATE=gsub("2020-08-26","2020-08-25", allData$DATE)
allData$PIG_DATE=gsub("2020-08-26","2020-08-25",allData$PIG_DATE)

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
#save(allData, file="data_clean/allData_minmeans_0_nomalize_T")

######
#MAYBE JUST USE META DIRECTLY, IS ALREADY MATCHED WITH FQ
######

#1) get metadata
#this is what we select FQs from
justMeta=allData[-6]
DF_meta=data.frame(justMeta)



forMimarks=data.frame(sample_name=DF_meta$PIG_DATE, 
                      sample_title="",
                      bioproject_accession="",
                      organism="unclassified sequences; metagenomes; organismal metagenomes; pig gut metagenome",
                      collection_date=DF_meta$DATE,
                      env_broad_scale="ENVO:01000998",
                      env_local_scale="ENVO:01001029",
                      env_medium="UBERON:0001988",
                      geo_loc_name="Denmark",
                      host="Sus scrofa domesticus",
                      lat_lon="not applicable",
                      OUA=DF_meta$OUA,
                      SOW=DF_meta$SOW,
                      group=DF_meta$BGROUP,
                      time=DF_meta$Week,
                      dummy=runif(NROW(DF_meta)))


write.table(file = "forMimarks.tsv",x = forMimarks,row.names = F,sep="\t")



#2) get FQs

#find fq parent path
FQ_path="C:/Users/milst/OneDrive - Danmarks Tekniske Universitet/Skrivebord/raw_data_GUDP_all/raw_data_GUDP_all/GUDP_opdeling/"

#get both subfolders
first="seq55"
second="seq56"

#read both
FIRST0  = dir(path = paste(FQ_path,first,sep=""))
SECOND0 = dir(path = paste(FQ_path,second,sep=""))

#combine both runs
BOTH_FQ=c(FIRST0,SECOND0)

#remove illumina endings and get only samples
FQ_SAMPLES=((gsub("_R._001.fastq.gz","", BOTH_FQ)))

#3) fish out sample name from meta, but fix
#starts with X
#has . instead of - in date
DF_meta$FULL_FQ=gsub("\\.","-",gsub("^X","",DF_meta$FULL))

DF_meta$R1=""
DF_meta$R2=""

i=1
for( i in 1:NROW(DF_meta)) {
  FQ_matches=grep(DF_meta$FULL_FQ[i],  BOTH_FQ, value = T) 
  if(length(FQ_matches)!=2) {
    print(DF_meta$FULL_FQ[i])

  } else {
    fwd=grep("R1",FQ_matches, value = T)
    rev=grep("R2",FQ_matches, value = T)
    DF_meta$R1[i]=fwd
    DF_meta$R2[i]=rev
  }
}

forSRA=data.frame(ID=DF_meta$PIG_DATE,
                  filename=DF_meta$R1,
                  filename2=DF_meta$R2)

write.table(x = forSRA,file = "forSRA.csv", row.names = F,sep=";")

#data.frame(head(DF_meta$FULL), head(FQ_SAMPLES))
