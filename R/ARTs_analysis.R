# Analyzing results from ARTs model
## Written by Sarah P. Flangan (spflanagan.phd@gmail.com)
### Date Started: 27 April 2017

setwd("~/Projects/ARTs/results/")
library(RColorBrewer)

## ---- CourterConditionalData
cc.sum<-read.delim("courter-conditional_summary.txt")
cc.traits<-read.delim("courter-conditional_traits.txt")
cc.popdyn<-read.delim("courter-conditional_popdyn.txt")

## ---- end

## ---- RandomMating
rm.sum<-read.table("random_mating_summary.txt",header=T)
plot(c(1,nrow(rm.sum)),c(0,1),type='n')
for(i in 7:ncol(rm.sum)){ points(rownames(rm.sum),rm.sum[,i],type="l")}

plot(rm.sum$CourterFreq, type="l")

rm.traits<-read.delim("random_mating_traits.txt")

## ---- end

## ---- CourterTrait
ct.sum<-read.table("courter_summary.txt",header=T)
ct.traits<-read.delim("courter_traits.txt")

plot(ct.sum$CourterFreq,type="l") #courters go to fixation
#do courter alleles go to fixation?
ct.qtl<-as.character(scan("courter_qtlinfo.txt",sep='\t',what=character()))
pref.qtl<-ct.qtl[grep("Pref",ct.qtl)]
crtr.qtl<-ct.qtl[grep("Courter",ct.qtl)]
#make them compatible with marker loci
pref.qtl<-gsub("PrefQTL(\\d.\\d+)","Marker\\1",pref.qtl)
crtr.qtl<-gsub("CourterQTL(\\d.\\d+)","Marker\\1",crtr.qtl)
pref.qtl.af<-ct.sum[,pref.qtl]
crtr.qtl.af<-ct.sum[,crtr.qtl]
marker.af<-ct.sum[,!(colnames(ct.sum) %in% pref.qtl) & !(colnames(ct.sum) %in% crtr.qtl)]
nmarkers<-ncol(marker.af)-6
bgcol<-gray.colors(nmarkers)
bgcol<-sample(bgcol)
ctcolsc<-colorRampPalette(c("white","navy blue"))
ctcol<-ctcolsc(ncol(crtr.qtl.af))
pfcolsc<-colorRampPalette(c("light pink","dark violet"))
pfcol<-pfcolsc(ncol(pref.qtl.af))

#plot
plot(c(1,nrow(ct.sum)),c(0,1),type='n')
for(i in 7:ncol(marker.af)){ points(rownames(marker.af),marker.af[,i],type="l",col=bgcol[i])}
for(i in 7:ncol(pref.qtl.af)){ points(rownames(pref.qtl.af),pref.qtl.af[,i],type="l",col=pfcol[i])}
for(i in 7:ncol(crtr.qtl.af)){ points(rownames(crtr.qtl.af),crtr.qtl.af[,i],type="l",col=ctcol[i])}

hist(unlist(crtr.qtl.af[nrow(crtr.qtl.af),]))
hist(unlist(pref.qtl.af[nrow(pref.qtl.af),]))
## ---- end

## ---- Courters4QTL
ct.sum<-read.table("courter_q4_summary.txt",header=T)
ct.traits<-read.delim("courter_q4_traits.txt")

plot(ct.sum$CourterFreq,type="l") #courters go to fixation
#do courter alleles go to fixation?
ct.qtl<-read.delim("courter_q4_qtlinfo.txt",sep='\t')
pref.qtls<-paste("Marker",ct.qtl[grep("Pref",colnames(ct.qtl))],sep="")
pref.qtl.af<-ct.sum[,pref.qtls]
crtr.qtls<-paste("Marker",ct.qtl[grep("Courter",colnames(ct.qtl))],sep="")
crtr.qtl.af<-ct.sum[,crtr.qtls]
#make them compatible with marker loci
marker.af<-ct.sum[,grep("Marker",colnames(ct.sum))]
nmarkers<-ncol(marker.af)
bgcol<-gray.colors(nmarkers)
bgcol<-sample(bgcol)
ctcolsc<-colorRampPalette(c("white","navy blue"))
ctcol<-ctcolsc(ncol(crtr.qtl.af))
pfcolsc<-colorRampPalette(c("light pink","dark violet"))
pfcol<-pfcolsc(ncol(pref.qtl.af))

#plot
plot(c(1,nrow(ct.sum)),c(0,1),type='n')
for(i in 1:ncol(marker.af)){ points(rownames(marker.af),marker.af[,i],type="l",col=bgcol[i])}
for(i in 1:ncol(pref.qtl.af)){ points(rownames(pref.qtl.af),pref.qtl.af[,i],type="l",col=pfcol[i])}
for(i in 1:ncol(crtr.qtl.af)){ points(rownames(crtr.qtl.af),crtr.qtl.af[,i],type="l",col=ctcol[i])}

hist(unlist(crtr.qtl.af[nrow(crtr.qtl.af),]))
hist(unlist(pref.qtl.af[nrow(pref.qtl.af),]))
## ---- end

## ---- ParentTrait
pt.sum<-read.table("parent_summary.txt",header=T)
pt.traits<-read.delim("parent_traits.txt")

plot(pt.sum$ParentFreq,type="l") 
#do parent alleles go to fixation?
pt.qtl<-read.delim("parent_qtlinfo.txt",sep='\t')
#make them compatible with marker loci
parent.qtls<-paste("Marker",pt.qtl[grep("Parent",colnames(pt.qtl))])
parent.qtl.af<-pt.sum[,parent.qtls]
marker.af<-pt.sum[,!(colnames(pt.sum) %in% parent)]
nmarkers<-ncol(marker.af)-6
bgcol<-gray.colors(nmarkers)
bgcol<-sample(bgcol)
ptcolsc<-colorRampPalette(c("tan","forest green"))
ptcol<-ptcolsc(ncol(parent.qtl.af))

#plot
plot(c(1,nrow(pt.sum)),c(0,1),type='n')
for(i in 7:ncol(marker.af)){ points(rownames(marker.af),marker.af[,i],type="l",col=bgcol[i])}
for(i in 7:ncol(parent.qtl.af)){ points(rownames(parent.qtl.af),parent.qtl.af[,i],type="l",col=ptcol[i])}

hist(unlist(parent.qtl.af[nrow(parent.qtl.af),]))
## ---- end

setwd("~/Projects/ARTs/results/")
## ---- FrequencyDependentPreference
pt.sum<-read.table("fdcourter_summary.txt",header=T)
pt.traits<-read.delim("fdcourter_traits.txt")

plot(pt.sum$CourterFreq,type="l") 
#do parent alleles go to fixation?
ct.qtl<-read.delim("fdpref_qtlinfo.txt",sep='\t')
#make them compatible with marker loci
pref.qtl<-paste("Marker",ct.qtl[grep("Pref",colnames(ct.qtl))],sep="")
pref.qtl.af<-pt.sum[,pref.qtls]
crtr.qtl<-paste("Marker",ct.qtl[grep("Courter",colnames(ct.qtl))],sep="")
crtr.qtl.af<-pt.sum[,crtr.qtls]
pref.qtl<-gsub("PrefQTL(\\d.\\d+)","Marker\\1",pref.qtl)
crtr.qtl<-gsub("CourterQTL(\\d.\\d+)","Marker\\1",crtr.qtl)
pref.qtl.af<-pt.sum[,pref.qtl]
crtr.qtl.af<-pt.sum[,crtr.qtl]
marker.af<-pt.sum[,!(colnames(pt.sum) %in% pref.qtl) & !(colnames(pt.sum) %in% crtr.qtl)]
nmarkers<-ncol(marker.af)-6
bgcol<-gray.colors(nmarkers)
bgcol<-sample(bgcol)
ctcolsc<-colorRampPalette(c("white","navy blue"))
ctcol<-ctcolsc(ncol(crtr.qtl.af))
pfcolsc<-colorRampPalette(c("light pink","dark violet"))
pfcol<-pfcolsc(ncol(pref.qtl.af))

#plot
plot(c(1,nrow(pt.sum)),c(0,1),type='n')
for(i in 7:ncol(marker.af)){ points(rownames(marker.af),marker.af[,i],type="l",col=bgcol[i])}
for(i in 7:ncol(pref.qtl.af)){ points(rownames(pref.qtl.af),pref.qtl.af[,i],type="l",col=pfcol[i])}
for(i in 1:ncol(crtr.qtl.af)){ points(rownames(crtr.qtl.af),crtr.qtl.af[,i],type="l",col=ctcol[i])}

hist(unlist(crtr.qtl.af[nrow(crtr.qtl.af),]))
hist(unlist(pref.qtl.af[nrow(pref.qtl.af),]))
## ---- end