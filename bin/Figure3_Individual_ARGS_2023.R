#####
#Mikael Lenz Strube & Katrine Tams
#A script to make timelineplots of individual ARGs on each sampling date, colored by treatment
#####rm(list = ls())

rm(list = ls())
source("bin/functions.R")

#Plot tp png file
makePNG=F

#Remove specific genes
remove_gene=TRUE
genes_to_remove=c("TolC1", "VatE", "water") # TolC is a e.coli housekeeoing gene, not an ARG, and VatE had technical issues on one of the chips squewing week 10 

library(vegan)
library(RColorBrewer)
library(beeswarm)
library(dunn.test)

allData=CleanData_qPCR(PATH="data/fluidigm_runs/",removeNul = F,condenseB = T)

#Removing unvanted gene 
if(remove_gene){
  allData3=allData
  allData3$DATA=allData3$DATA[,-which(colnames(allData$DATA) %in% genes_to_remove)]
  allData=allData3
}

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

#Groups in choosen groups 
uniq_groups=levels(factor(unique(allData$BGROUP)))

#Colors for plotting
pal=c("red", "purple", "blue", "black")

names(pal)=c("A","B_T","B_U","C")


cols=pal[factor(allData$BGROUP)]

levels(factor(unique(allData$Week)))

geneWeekGroupDF=data.frame(GENEINDX="", GENE="",WEEK="",P="", LARGEST="")

geneIndx=1


runningCounter=1
for( geneIndx in 1:NCOL(allData$DATA) ) {
  
  geneDummy=data.frame(WEEK=allData$Week,GROUP=allData$BGROUP,VALS=allData$DATA[,geneIndx])
  
  wk="Week 00"
  for(wk in levels(factor(unique(allData$Week)))) {
    weekGeneDummy=subset(geneDummy, WEEK==wk)
    weekGeneMedians=aggregate(VALS~GROUP,data = weekGeneDummy, median)
    
    KW_geneWeek=kruskal.test(VALS~GROUP,weekGeneDummy)
    dummy_vec=c(geneIndx,
                colnames(allData$DATA)[geneIndx],
                wk,
                KW_geneWeek$p.value,
                weekGeneMedians[order(weekGeneMedians$VALS, decreasing = T),][1,1]
    )  
    geneWeekGroupDF=rbind(geneWeekGroupDF,dummy_vec)
    runningCounter=runningCounter+1
  }
  
}

geneWeekGroupDF$adjP=p.adjust(geneWeekGroupDF$P, method = "BH")

geneWeekGroupDF=geneWeekGroupDF[-1,]

geneWeekGroupDF$GENE=factor(geneWeekGroupDF$GENE)
geneWeekGroupDF$WEEK=factor(geneWeekGroupDF$WEEK)
geneWeekGroupDF$P=as.numeric(geneWeekGroupDF$P)
geneWeekGroupDF$adjP=as.numeric(geneWeekGroupDF$adjP)


sigGenes=c()
i = unique(geneWeekGroupDF$GENE)[9]
for(i in unique(geneWeekGroupDF$GENE)) {
  tempGene=subset(geneWeekGroupDF, GENE==i)
  
  if(length(which(ifelse(tempGene$adjP<0.05,1,0)==1))>1) {
    sigGenes=c(sigGenes,i)
    print(tempGene)
  }
}

#sigGenes=unique(geneWeekGroupDF[which(geneWeekGroupDF$adjP<0.05),2])


#####
#sigGenes
#####

#Open png generater
if(makePNG) png(filename = "Results/Fig_3_sing_ARGs_timeline_oneplot_longnames2023.png", width = 7.4, height = 8,units = "in" , res=600)


#setup plot parameters, e.g. a 3x3 matrix with very short last row
par(mar=c(0,2.5,1.5,0), mgp=c(1.4,0.4,0), font.lab=2)
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15),nrow = 5,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.5,0.5,0.5,0.5,0.15), widths = c(0.3,0.3,0.3))


plotCounter=1
for(sGene in sigGenes) {
  byGeneSum=aggregate(allData$DATA[,sGene]~allData$Week + allData$BGROUP, 
                      FUN = median); 
  
  colnames(byGeneSum)=c("WEEK","GROUP","MEDIANS")
  byGeneSum$WEEK=factor(byGeneSum$WEEK)
  
  
  YLIM=c(-max(byGeneSum$MEDIANS)/5,max(byGeneSum$MEDIANS)*1.5)
  if(plotCounter >7) {
    beeswarm(allData$DATA[,sGene]~allData$Week, corral="wrap", ,main=sGene,
           pwcol=cols, pch=16, cex=0.5, las=2,xlab="",ylab="", ylim=YLIM)
  } else {
    beeswarm(allData$DATA[,sGene]~allData$Week, corral="wrap", xaxt="n",main=sGene,
             pwcol=cols, pch=16, cex=0.5, las=2,xlab="",ylab="", ylim=YLIM)
  }
  
  lines(x = c(2,2),y = c(0,10), col="grey50", lty=2,lwd=4)
  for(grp in unique(byGeneSum$GROUP)) {
    lines(MEDIANS~WEEK,subset(byGeneSum, GROUP==grp), col=pal[grp],lwd=2)
  }
  
  geneForPlot=subset(geneWeekGroupDF, GENE==sGene & P<0.05)
  
  fullGeneDat=data.frame(VAL=allData$DATA[,sGene],WEEK=allData$Week, GROUP=allData$BGROUP)
  
  j= geneForPlot$WEEK[1]
  for(j in geneForPlot$WEEK) {
    sigWeekGene=subset(fullGeneDat, WEEK==j)
    #DUNN=dunn.test(x =sigWeekGene$VAL,g = sigWeekGene$GROUP,kw = T, table = F )
    KW=kruskal.test(x =sigWeekGene$VAL,g = sigWeekGene$GROUP)
    if(KW$p.value<0.05) {
      text(which(j==levels(geneForPlot$WEEK)),0.9*-max(byGeneSum$MEDIANS)/5,"*", cex=2)
    }
  }
  
  plotCounter=plotCounter+1
}

#Make empty plot 
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")

#make minimal plot and horizontal legend
par(mar=c(0,0,0,0))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
#Legend depending on groups
if(length(uniq_groups)==2){
  legend(x = "center",legend = c("Treated", "Untreated"), col = pal,pch = 16, horiz = F)
}else{    
  legend_groups=levels(factor(uniq_groups))
  legend(x = "bottomright",legend = c("Treated group", "Mix group treated","Mix group untreated", "Untreated", "Time of treatment", "Significant difference"), col = c(pal,"black","black"),pch = c(16, 16, 16, 16, NA, 42), lty = c(NA, NA, NA, NA, 2, NA),  horiz = F, cex=1.7)
}


if(makePNG) dev.off()
###


