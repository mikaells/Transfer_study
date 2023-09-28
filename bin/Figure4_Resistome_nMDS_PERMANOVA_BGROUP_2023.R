#####
#Mikael Lenz Strube & Katrine Tams
#A script to make NMDS plots for each sampling date of the 16S data colored by treatment time and a PERMANOVA analysis of the same groups

rm(list = ls())
source("bin/functions.R")

library(vegan)
library(beeswarm)
library(RColorBrewer)
library(phyloseq)

addEnvfit=F
makePNG=F

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
load("data_clean/PERMANOVA_RESISTOME_BGROUP_with_week_0")
#Or run the loop to create the PERMANOVA values
# date=uniq_dates[1]
# for (date in uniq_dates){
#   sample_date=allData$DATE[which(allData$DATE==date)[1]]
#   sample_week=allData$Week[which(allData$DATE==date)[1]]
# 
#   print(sample_date)
#   print(sample_week)
# 
#   #Index for samples at given date
#   dateIndex=which(allData[["DATE"]]==date)
# 
#   #To ensure there is at least two groups to run the PERMANOVA
#   #Creating dataframes for the given date
#   dateData=allData[["DATA"]][dateIndex,]
#   dateOUA=allData[["OUA"]][dateIndex]
#   dateBGROUP=allData[["BGROUP"]][dateIndex]
#   dateNote=allData[["Note"]][dateIndex]
# 
# 
#   #Calculating PERMANOVA for given date
#   ADO=adonis2(dateData~factor(dateBGROUP), parallel = 2, permutations = 200000) #
# 
#   #Calculating the betadispersion between the groups
#   DISP=anova(betadisper(d = vegdist(dateData),group = factor(dateBGROUP)))
#   #print(betadisper(d = vegdist(dateData),group = factor(dateBGROUP)))
# 
# 
#   Rsquared=c(Rsquared, signif(ADO$R2[1],4))
#   Rdates=c(Rdates, date)
#   Rweeks=c(Rweeks, sample_week)
#   pValue=c(pValue, signif(ADO$`Pr(>F)`[1]))
#   pdates=c(pdates, date)
#   pweeks=c(pweeks, sample_week)
# }
# 
# #pValues are adjusted to the Benjamini & Yekutieli (2001) ("BY") https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust
# pAdjusted=p.adjust(p=pValue, method = "BY")
# pvalue_condensed=formatC(pAdjusted, format = "e", digits = 1)
# 
# PERMANOVA_values=as.data.frame(cbind(pweeks, Date=pdates,Rsquared=Rsquared, pValue, pAdjusted, pvalue_condensed))
# save(PERMANOVA_values,file = paste("data_clean/PERMANOVA_RESISTOME_BGROUP_with_week_0"))

#write.table(PERMANOVA_values, row.names = F, col.names = T, file = paste("Results/Resistome_PERMANOVA_values.txt"))




#open plot generators
if(makePNG) png(filename = paste("Results/Fig4_Resistome_MainFig_spider_BGROUP_week_0.png"), width = 7.4, height = 6,units = "in" , res=600)


#setup plot parameters, e.g. a 3x3 matrix with very short last row
par(mar=c(2.5,2.5,0.5,.5), mgp=c(1.4,0.4,0), font.lab=2)
m <- matrix(c(1,2,3,4,5,6,7,8,9, 10, 11, 12),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(0.6,0.65,0.65,0.2), widths = c(0.65,0.6,0.6,0.6))

#Creating colors for nMDS plot
pal=c("red", "purple", "blue", "black")
cols=pal[factor(allData$BGROUP)]
pchs_treatment=c(1,16)
pchs=pchs_treatment[factor(allData$Note)]

plot_count=0
#date=uniq_dates[1]
for (date in uniq_dates){
  sample_date=allData$DATE[which(allData$DATE==date)[1]]
  sample_week=allData$Week[which(allData$DATE==date)[1]]
  print(sample_date)
  print(sample_week)
  
  PERM_dateIndex=which(PERMANOVA_values$Date==date)
  dateIndex=which(allData[["DATE"]]==date)
  
  #Date for given date
  dateData=allData[["DATA"]][dateIndex,]
  dateOUA=allData[["OUA"]][dateIndex]
  dateBGROUP=allData[["BGROUP"]][dateIndex]
  
  #Calculating nMDS
  NMDS=metaMDS(dateData, parallel=2, trace = 0)
  NMDSscores=scores(NMDS, display = "sites")
  # distFromCenter = (sqrt(NMDSscores[, 1]^2 + NMDSscores[, 2]^2))
  # #plot(distFromCenter)
  # OUTS = as.numeric((which(distFromCenter > 5 *
  #                            mean(distFromCenter))))
  #distFromCenter = (sqrt(NMDSscores[, 1]^2 + NMDSscores[, 2]^2))
  distFromCenter1 = (sqrt(NMDSscores[, 1]^2))
  distFromCenter2 = (sqrt(NMDSscores[, 2]^2))
  #plot(distFromCenter)
  OUTS = as.numeric((which(distFromCenter1 > 5 *
                             mean(distFromCenter1))))
  OUTS = c(OUTS, as.numeric((which(distFromCenter2 > 5 *
                                     mean(distFromCenter2)))))
  OUTS=unique(OUTS)
  
  if (length(OUTS) > 0) {
    cat(paste("Sample(s)\n"))
    cat(paste("\t", OUTS))
    cat("\nis/are suspecious, consider removing them with 'SampleRemover'.\n")
    
    cat(paste("Sample(s) removed\n"))
    dateData=dateData[-c(OUTS),]
    dateBGROUP=dateBGROUP[-c(OUTS)]
    dateOUA=dateOUA[-c(OUTS)]
    dateIndex=dateIndex[-c(OUTS)]
    # 
    NMDS=metaMDS(dateData, parallel=1, trace = 0)
    
    NMDSscores = scores(NMDS, display = "sites")
  }
  
  #Values for plot title
  pvalue_for_title=as.numeric(PERMANOVA_values$pvalue_condensed[PERM_dateIndex])
  Rsquared_titel=as.numeric(PERMANOVA_values$Rsquared[PERM_dateIndex])
  
  datepchs=pchs[dateIndex]
  
  plot_count=plot_count+1
  par(mar=c(0.1,0.1,0.1,0.1))
  plot(scores(NMDS, display = "sites"), col="white", pch=16, 
       ylim=c(min(scores(NMDS, display="sites")[,2]),1.3*max(scores(NMDS, display="sites")[,2])), xaxt='n',yaxt='n', xlab="",ylab="")    
  abline(v=0, col="grey")
  abline(h=0, col="grey")
  points(scores(NMDS, display="sites"), col=cols[dateIndex], pch=datepchs, cex = 1)
  #text(scores(NMDS, display="sites"),allData$PIG_DATE[dateIndex])
  #text(scores(NMDS, display="sites"))
  
  
  #text(scores(NMDS, display="sites"),labels = dateBGROUP)
  ordiellipse(NMDS,groups =factor(dateBGROUP),col = pal, kind="se", conf=0.95)
  ordispider(NMDS,groups =factor(dateBGROUP),col = pal)
  #ordiellipse(NMDS,groups =factor(dateOUA),col = c("red", "blue") , kind="se", conf=0.95 , lty = 2)
  
  
  legend(x = "top",legend = bquote(.(sample_week) ~"|"~ R^2~"="~ .(Rsquared_titel*100) ~ "% |"~" P=" ~ .(pvalue_for_title)), box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==1)legend(x = "topleft", legend ="a)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==2)legend(x = "topleft", legend ="b)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==3)legend(x = "topleft", legend ="c)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==4)legend(x = "topleft", legend ="d)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==5)legend(x = "topleft", legend ="e)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==6)legend(x = "topleft", legend ="f)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==7)legend(x = "topleft", legend ="g)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==8)legend(x = "topleft", legend ="h)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  if(plot_count==9)legend(x = "topleft", legend ="i)", box.lwd = 0,box.col = "transparent",bg = "transparent")
  
  if(addEnvfit) {
    ENV=envfit(NMDS,env = dateData)
    top10minR=min(sort(ENV$vectors$r, decreasing = T)[1:10])
    sigEnvIndx=which(ENV$vectors$r>top10minR)
    ENV$vectors$arrows=ENV$vectors$arrows[sigEnvIndx,]
    ENV$vectors$pvals=ENV$vectors$pvals[sigEnvIndx]
    ENV$vectors$r=ENV$vectors$r[sigEnvIndx]
    plot(ENV)
  }
  
}

#make minimal plot and legend
par(mar=c(0.1,0.1,0.1,0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
box(col="transparent")

par(mar=c(0.1,0.1,0.1,0.1))
plot(1, type = "n", axes=FALSE, xlab="", ylab="")
box(col="transparent")
#legend(x = "top",   pch=16, cex=.9, horiz = T, text.width=c(.03, 0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015), legend = c("Weeks", "3", "4",  "5", "6","7", "10", "13", "16"), col = c("transparent",pal_Date),box.col = "transparent",bg = "transparent")
legend(x = "center",  pch=16 ,text.width=c(.05,.09,.09,.09), cex=.9, horiz = T, legend = c( "A", "B_T", "B_U", "C"), col = c(pal),box.col = "transparent",bg = "transparent")



#close plot generators
if(makePNG) dev.off()

#For SIF 
# ADO_all_time_OUA=adonis2(allData$DATA~factor(allData$OUA)*factor(allData$SAMPLEWEEK), parallel = 18, permutations = 200000)
# save(ADO_all_time_OUA,file = "data_clean/PERMANOVA_adonis2_alltime_OUA")
# ADO_all_time_treattime=adonis2(allData$DATA~factor(allData$WEEK)*factor(allData$SAMPLEWEEK), parallel = 18, permutations = 200000)
# save(ADO_all_time_OUA,file = "data_clean/PERMANOVA_adonis2_alltime_treattime")


