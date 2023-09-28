#####
#Mikael Lenz Strube & Katrine Tams
#A script to make timelineplots of individual ARGs on each sampling date, colored by treatment
#####rm(list = ls())

rm(list = ls())
source("bin/functions.R")

#Plot tp png file
makePNG=F
calcSig=F
addSig=T

#Remove specific genes
remove_gene=TRUE
genes_to_remove=c("TolC1", "VatE", "water", "X16s", "X16s.1", "X16s.2") # TolC is a e.coli housekeeoing gene, not an ARG, and VatE had technical issues on one of the chips squewing week 10 

library(vegan)
library(RColorBrewer)
library(beeswarm)
library(dunn.test)
library(naturalsort)

allData=CleanData_qPCR(PATH="data/fluidigm_runs/",removeNul = F,condenseB = T)

rowSums(allData$DATA)

#Removing unvanted gene 
if(remove_gene){
  allData3=allData
  allData3$DATA=allData3$DATA[,-which(colnames(allData$DATA) %in% genes_to_remove)]
  allData=allData3
}

rowSums(allData$DATA)

#Alter sampling date for plotting purpose 
allData$DATE=gsub("2020-08-26","2020-08-25", allData$DATE)
allData$PIG_DATE=gsub("2020-08-26","2020-08-25",allData$PIG_DATE)

##Adding treatment data 
#Reading in meta data 
meta_data <- read.delim("data/meta_data.txt")
meta_data$DATE=as.Date(meta_data$DATE,format = "%d-%m-%Y" )
meta_data$pig_date=str_c(meta_data$PIG, '_', meta_data$DATE)

#Creating placeholder for metadata
allData$OUA=allData$PIG
allData$Stald=allData$PIG
allData$TreatTime=allData$PIG
allData$Week=allData$PIG
allData$BGROUP=allData$PIG

#Looping over data to insert metadata
for(i in 1:length(allData$PIG_DATE)){
  if(isTRUE(allData$PIG_DATE[i] %in% meta_data$pig_date)){
    metaIndex=which(meta_data$pig_date==allData$PIG_DATE[i])
    allData$OUA[i]=meta_data$OUA[metaIndex]
    allData$Stald[i]=meta_data$Stald[metaIndex]
    allData$TreatTime[i]=meta_data$Behandlingsdato[metaIndex]
    allData$Week[i]=meta_data$Week[metaIndex]
    allData$BGROUP[i]=meta_data$GROUP_B[metaIndex]
    
  }else{
    allData$OUA[i]=NA
    allData$Stald[i]=NA
    allData$TreatTime[i]=NA
    allData$Week[i]=NA
    allData$BGROUP[i]=NA
  }
} 



#Importing names for plots and replacing colnames 
Gene_names <- read_excel("data/Gene_names.xlsx", 
                         range = "A1:C86", col_names = T)
for(i in 1:NCOL(allData$DATA)){
  if(colnames(allData$DATA)[i] %in% Gene_names[["Names_on_chip"]]){
    index=which(Gene_names[["Names_on_chip"]]==colnames(allData$DATA)[i])
    colnames(allData$DATA)[i]=Gene_names[["Type_name(primer)"]][index]    
  }
  
}



#Calculating the sum of each gene for sumplot 
resistance_sums=rowSums(allData$DATA)


#####
#OUA comparisons
#####

geneDummyOUA_unsorted=data.frame(WEEK=allData$Week,GROUP=allData$OUA,VALS=resistance_sums)
geneDummyOUA=geneDummyOUA_unsorted[naturalorder(geneDummyOUA_unsorted$WEEK),]


WeekOUADF=data.frame(WEEK="",P="", LARGEST="")

wk="Week 00"
runningCounter=0
for(wk in unique(geneDummyOUA$WEEK)) {
  weekOUADummy=subset(geneDummyOUA, WEEK==wk)
  weekOUAMedians=aggregate(VALS~GROUP,data = weekOUADummy, median)
  
  if(wk=="Week 00") {
    
    dummy_vec=c(wk,  1,
                weekOUAMedians[order(weekOUAMedians$VALS, decreasing = T),][1,1]
    ) 
  } else {
    KW_OUAWeek=kruskal.test(VALS~GROUP,weekOUADummy)
    dummy_vec=c(wk,  KW_OUAWeek$p.value,
                weekOUAMedians[order(weekOUAMedians$VALS, decreasing = T),][1,1]
    )  
  }
  WeekOUADF=rbind(WeekOUADF,dummy_vec)
  runningCounter=runningCounter+1
  
}

WeekOUADF$adjP=p.adjust(WeekOUADF$P, method = "BH")

WeekOUADF=WeekOUADF[-1,]

WeekOUADF$WEEK=factor(WeekOUADF$WEEK)
WeekOUADF$P=as.numeric(WeekOUADF$P)
WeekOUADF$adjP=as.numeric(WeekOUADF$adjP)

#####
#BGROUP comparisons
#####


geneDummyBGR_unsorted=data.frame(WEEK=allData$Week,GROUP=allData$BGROUP,VALS=resistance_sums)
geneDummyBGR=geneDummyBGR_unsorted[naturalorder(geneDummyBGR_unsorted$WEEK),]


WeekBGRDF=data.frame(WEEK="",P="", LARGEST="",
                     "AvsB_T"="",  "AvsB_U" ="", "B_TvsB_U"="", "AvsC"="", "B_TvsC"="",   "B_UvsC"="" )

wk="Week 03"
runningCounter=0
for(wk in unique(geneDummyBGR$WEEK)) {
  weekBGRDummy=subset(geneDummyBGR, WEEK==wk)
  weekBGRMedians=aggregate(VALS~GROUP,data = weekBGRDummy, median)
  
  if(wk=="Week 00") {
    
    dummy_vec=c(wk,  1,
                weekBGRMedians[order(weekBGRMedians$VALS, decreasing = T),][1,1],
                rep(1,6)
    ) 
  } else {
    KW_BGRWeek=kruskal.test(VALS~GROUP,weekBGRDummy)
    DUNN=dunn.test(x = weekBGRDummy$VALS, g = weekBGRDummy$GROUP, table = F)
    dummy_vec=c(wk,  KW_BGRWeek$p.value,
                weekBGRMedians[order(weekBGRMedians$VALS, decreasing = T),][1,1],
                DUNN$P
    )  
  }
  WeekBGRDF=rbind(WeekBGRDF,dummy_vec)
  runningCounter=runningCounter+1
  
}

WeekBGRDF$adjP=p.adjust(WeekBGRDF$P, method = "BH")

WeekBGRDF=WeekBGRDF[-1,]

WeekBGRDF$WEEK=factor(WeekBGRDF$WEEK)
WeekBGRDF$P=as.numeric(WeekBGRDF$P)
WeekBGRDF$adjP=as.numeric(WeekOUADF$adjP)

#####
#plot OUA
#####

#Open png generater
if(makePNG) png(filename = "Results/Fig2_medSums_2023.png", width = 7.4, height = 6,units = "in" , res=600)


par(mar=c(5,3,1.5,0), mgp=c(1.6,0.5,0), font.lab=2, mfrow=c(1,2))
beeswarm(resistance_sums ~ allData$Week, pch=16, ylim=c(0,1.3), 
         cex=0.5, pwcol = ifelse(as.numeric(allData$OUA)==1,1,2), las=2, 
         xlab="", ylab="Median sum of genes per 16S gene")
abline(v = 1.5)

legend(x = "topright",legend = c("Treated", "Untreated"), col=c(2,1),pch=16)
OUA_plot=0
for(OUA_plot in unique(geneDummyOUA$GROUP)) {
  dummyOUAplot=subset(geneDummyOUA, GROUP==OUA_plot)
  OUA_week=aggregate(VALS~WEEK, dummyOUAplot, FUN=median)
  
  #shift line one week if treated
  #no treated in week0
  adjuster=ifelse(OUA_plot==0,1,0)
  
  lines(as.numeric(factor(OUA_week$WEEK))+adjuster, OUA_week$VALS, lwd=3, 
        col = ifelse(as.numeric(unique(dummyOUAplot$GROUP))==1,1,2))
}

if(addSig) {
  text(c(2:6),rep(1.2,4),rep("*",4), cex=2)
}

if(calcSig) {
  for(j in 1:NROW(WeekOUADF)) {
    if(WeekOUADF$adjP[j] <0.05 ) {
      text(WeekOUADF[j,1],1.2,"*",cex=3)
    } 
  }
}
#####
#plot GRPB
#####

#Colors for plotting
pal=c("red", "purple", "blue", "black")
names(pal)=c("A","B_T","B_U","C")

pal[allData$BGROUP]

#par(mar=c(3,2.5,1.5,0), mgp=c(1.4,0.4,0), font.lab=2)
beeswarm(resistance_sums ~ allData$Week, pch=16, ylim=c(0,1.3), 
         cex=0.5, pwcol = pal[allData$BGROUP], las=2, 
         xlab="", ylab="Median sum of genes per 16S gene")
abline(v = 1.5)

for(BGR_plot in unique(geneDummyBGR$GROUP)) {
  dummyBGRplot=subset(geneDummyBGR, GROUP==BGR_plot)
  BGR_week=aggregate(VALS~WEEK, dummyBGRplot, FUN=median)
  
  
  lines(as.numeric(factor(BGR_week$WEEK)), BGR_week$VALS, lwd=3, 
        col = pal[unique(dummyBGRplot$GROUP)])
  legend(x = "topright",
         legend = c("A", "B_T","B_U", "C"), col = pal,pch = c(16, 16, 16, 16))
  
}

if(addSig) {
  text(c(2:7),rep(1.2,5),c("a","a","b","c","b","d"))
}

if(calcSig) {
  
  for(j in 1:NROW(WeekOUADF)) {
    if(WeekOUADF$adjP[j] <0.05 ) {
      text(WeekOUADF[j,1],1.2,"*",cex=3)
    } 
  }
}
if(makePNG) dev.off()


