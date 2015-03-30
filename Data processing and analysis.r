library (limma)

# read Targets and raw data
Targets<- readTargets("targets.txt",row.names="Names")
RG<-read.maimages(Targets$Filename, source="agilent")

# positive and negative control probes are removed
non_controls <- which(RG$genes$ControlType == 0)
dataset_no_controls<- RG[-non_controls,]

# Low expressed genes (below the 20th percentile) are removed
fake_MA<- MA.RG(dataset_no_controls)
try<- apply(fake_MA$A, 1, median, na.rm=T)
quantile(fake_MA$A,na.rm=T)
quantile(fake_MA$A,.20,na.rm=T)
jpeg("boxplot_filter.jpeg",width=6000, height=3000,quality=100,res=300)
boxplot(fake_MA$A)
abline(h=4.254393, col="red")
dev.off()
low_expressed<- which(try <= 4.254393)
final_dataset<-dataset_no_controls[-low_expressed,]

#Data are Loess and then quantile normalized
MA<-normalizeWithinArrays(final_dataset, method="loess")
MA.q<-normalizeBetweenArrays(MA, method="quantile")

#MA plots
for (i in 1:106){
	jpeg(paste("MA_plot_",i,".jpeg"),width=6000, height=3000,quality=100,res=300)
	par(mfrow=c(1,3))
   	plotMA(RG[,i], ylim=c(-6,6),status=RG$genes$ControlType)
	abline(h=0, col="red")
	plotMA(MA[,i], ylim=c(-6,6))
	abline(h=0, col="red")
	plotMA(MA.q[,i], ylim=c(-6,6))
	abline(h=0, col="red")
	dev.off()
}


# quality control process
library(arrayQualityMetrics)
arrayQualityMetrics(RG, outdir= "Quality_Assessment_Raw_data", do.logtransform=T)
arrayQualityMetrics(MA, outdir= "Quality_Assessment_Loess")
arrayQualityMetrics(MA.q, outdir= "Quality_Assessment_Quantile")

#density plots
jpeg("Density_plot.jpeg",width=6000, height=3000,quality=100,res=300)
par(mfrow=c(1,3))
plotDensities(RG)
plotDensities(MA)
plotDensities(MA.q)
dev.off()

#duplicate probes are summarized by average
summarize.malist <- function(DATA,by=c("max","avg"),na.rm=T,...){
	ratios <- DATA$M
	intensities <- DATA$A
	annot <- DATA$genes
	sorted <- tapply(1:nrow(ratios),annot$ProbeName,function(x) x)
	test <- lapply(sorted,function(x){
		res <- switch(by,
			max = {
				sel <- x[which.max(rowMeans(intensities[x,,drop=F],na.rm=na.rm))]
				return(list(ratios=ratios[sel,],intensities=intensities[sel,],annot=annot[sel,]))
			},
			avg = {
				return(list(ratios=colMeans(ratios[x,,drop=F]),intensities=colMeans(intensities[x,,drop=F]),annot=annot[x[1],,drop=F]))
			})
		return(res)
	})
	ratios.n <- t(sapply(test,function(x) x$ratios))
	intensities.n <- t(sapply(test,function(x) x$intensities))
	annot.n <- t(sapply(test,function(x) as.character(x$annot)))
	colnames(annot.n) <- colnames(annot)
	annot.n <- data.frame(annot.n)
	DATA2 <- DATA
	DATA2$M <- ratios.n
	DATA2$A <- intensities.n
	DATA2$genes <- annot.n
	return(DATA2)
}

dataset_summ <- summarize.malist(MA.q,by="avg")

#NA values are removed
M_values_from_final_dataset <- dataset_summ$M
NA_<- apply(M_values_from_final_dataset,1,function(x) length(which(is.na(x)))/length(x))
to_remove<-which(NA_ > 0)
dataset_ready<-dataset_summ[-to_remove,]


#annotation file is read in and a gene name is assigned to each probes. NB: We kept those genes with a unique gene symbol assigned (column 4 == 1)
annotation<-read.delim("annotation.txt",header=T)
annotation.s <- by(annotation,annotation$ProbeID,function(x) x)
annotation.s <- do.call(rbind,annotation.s)
dataset_complete_annotated.d <- dataset_ready[intersect(rownames(dataset_ready),rownames(annotation.s)),]
annotation.d <- annotation.s[rownames(dataset_complete_annotated.d),]

#duplicated genes are summarized keeping the one with the highest averaged intensity
summarize <- function(DATA,INDICES,by=c("max","avg"),return.unlisted=T,na.rm=T,intensities=NULL){	
	if(!is.null(intensities)){
		res <- switch(by,
			max = tapply(1:nrow(DATA),INDICES,function(x) return(as.numeric(x[which.max(rowMeans(intensities[x,,drop=F],na.rm=na.rm))]))),
			avg = by(DATA,INDICES,function(x) return(colMeans(x,na.rm=na.rm)))
		)
		res <- t(sapply(res,function(x) x))	
	}else{
		res <- switch(by,
			max = by(DATA,INDICES,function(x) return(x[which.max(rowMeans(x,na.rm=na.rm)),])),
			avg = by(DATA,INDICES,function(x) return(colMeans(x,na.rm=na.rm)))
		)
		res <- t(sapply(res,function(x) x))
		if(return.unlisted)
			res <- delist(res)
	}
	return(res)
}

genes_wanted<-summarize(dataset_complete_annotated.d$M,annotation.d[,2],by="max", intensities=dataset_complete_annotated.d$A)
genes_wanted2<-t(genes_wanted)
dataset_complete_annotated_final<-dataset_complete_annotated.d[genes_wanted2[,1],]

#A new target file is created which is needed for channel separation
Targets2<-Targets
names(Targets2)<-c("Hybe","Cy3","Cy5","FileName")
cambio1<-gsub(" ","",Targets2$Cy3)
cambio2<-gsub(" ","",Targets2$Cy5)
Targets2$Cy3<-cambio1
Targets2$Cy5<-cambio2

#Channel are separated
separated_channel<-RG.MA(dataset_complete_annotated_final)
red_intensities<-separated_channel$R
green_intensities<-separated_channel$G
red_intensities_log<-log2(red_intensities)
green_intensities_log<-log2(green_intensities)
colnames(green_intensities_log)<- Targets2$Cy3
colnames(red_intensities_log)<- Targets2$Cy5
dataset_separated_channel<- cbind(green_intensities_log,red_intensities_log)
rownames(dataset_separated_channel)<-dataset_complete_annotated_final[[2]]$ProbeName
dataset_separated_channel.d<-na.omit(dataset_separated_channel)

#we made a file where in the first column there are all the chemicals samples and in the second column their relatives control
samples<-read.delim("samples.txt",header=T)
samples.l <- list()
samples.l[[1]] <- lapply(samples[,1],function(x) which(colnames(dataset_separated_channel.d) == x))
names(samples.l[[1]]) <- samples[,1]
samples.l[[2]] <- lapply(samples[,2],function(x) which(colnames(dataset_separated_channel.d) == x))
names(samples.l[[2]]) <- samples[,2]

#Since there is no an equal number of replicates for each exposure group, controls for each exposure group are averaged
control.avg <- lapply(samples.l[[2]],function(x) rowMeans(dataset_separated_channel.d[,x]))
#data ratios is then calculated
ratios <- lapply(1:length(samples.l[[1]]),function(x) dataset_separated_channel.d[,samples.l[[1]][[x]]]-control.avg[[x]])
data.ratios <- do.call(cbind,ratios)

#dataset are written into a file
dataset_separated_channel_final<-cbind(rownames(genes_wanted2),dataset_separated_channel.d)
dataset_ratios_final<-cbind(rownames(genes_wanted2),data.ratios)
colnames(dataset_separated_channel_final)[[1]]<-"GeneName"
colnames(dataset_ratios_final)[[1]]<-"GeneName"
write.table(dataset_separated_channel_final,"Dataset_separated_channel.txt",sep="\t",quote=F,col.names=NA)
write.table(dataset_ratios_final,"Dataset_ratios.txt",sep="\t",quote=F,col.names=NA)

__________________________________________________________________________________________________________________________

# Single chemical loop analysis

# Design matrix is created and then a linear model is fitted (we chose a general "ref" but the real contrast matrix is calculated later)
design<- modelMatrix(Targets2, ref="A/ECtrl")
fit<- lmFit(dataset_complete_annotated_final, design)

#Here we make contrasts
contrast_matrix<-read.delim("contrast_matrix.txt",header=T,row.names=1)
fit_ChemVsChem<- contrasts.fit(fit,contrast_matrix)
fit2_ChemVsChem<-eBayes(fit_ChemVsChem)

# DE genes are retrieved for three different thresholds
comparison1_ChemVsChem<-c()
for (i in 1: length(colnames(contrast_matrix))){
	comparison1_ChemVsChem[[i]]<-topTable(fit2_ChemVsChem, p.value=0.01, number= 50000, genelist=fit2_ChemVsChem$genes, coef= colnames(contrast_matrix)[[i]], resort.by="t")
	
}

comparison5_ChemVsChem<-c()
for (i in 1: length(colnames(contrast_matrix))){
	comparison5_ChemVsChem[[i]]<-topTable(fit2_ChemVsChem, p.value=0.05, number= 50000, genelist=fit2_ChemVsChem$genes, coef= colnames(contrast_matrix)[[i]], resort.by="t")
	
}

comparison10_ChemVsChem<-c()
for (i in 1: length(colnames(contrast_matrix))){
	comparison10_ChemVsChem[[i]]<-topTable(fit2_ChemVsChem, p.value=0.1, number= 50000, genelist=fit2_ChemVsChem$genes, coef= colnames(contrast_matrix)[[i]], resort.by="t")
	
}

#We get genes up and down regulated and we retrieve their own annotation
comparison1_ChemVsChem.s <- lapply(comparison1_ChemVsChem,function(x) try(by(x,sign(x[,"logFC"]),function(y) y)))
comparison1_ChemVsChem_Genes<-lapply(comparison1_ChemVsChem.s,lapply,function(x){
	if(length(grep("^Error in",x[1])) == 0){
		ann <- annotation.d[match(x$ProbeName,annotation.d$ProbeID),]
	}
	else{
		ann <- NA
	}
	return(list(annotation=ann))
})

comparison5_ChemVsChem.s <- lapply(comparison5_ChemVsChem,function(x) try(by(x,sign(x[,"logFC"]),function(y) y)))
comparison5_ChemVsChem_Genes<-lapply(comparison5_ChemVsChem.s,lapply,function(x){
	if(length(grep("^Error in",x[1])) == 0){
		ann <- annotation.d[match(x$ProbeName,annotation.d$ProbeID),]
	}
	else{
		ann <- NA
	}
	return(list(annotation=ann))
})

comparison10_ChemVsChem.s <- lapply(comparison10_ChemVsChem,function(x) try(by(x,sign(x[,"logFC"]),function(y) y)))
comparison10_ChemVsChem_Genes<-lapply(comparison10_ChemVsChem.s,lapply,function(x){
	if(length(grep("^Error in",x[1])) == 0){
		ann <- annotation.d[match(x$ProbeName,annotation.d$ProbeID),]
	}
	else{
		ann <- NA
	}
	return(list(annotation=ann))
})




#here we get gene lists
com.up <- function(x){
	x <- lapply(x,do.call,what=cbind)
	up<-c()
	lo<-c()
	if(any(names(x)=="-1")){
		lo<-rbind(lo,as.matrix(cbind(x[["-1"]],"down")))
	}
	if(any(names(x)=="1")){
		up<-rbind(up,as.matrix(cbind(x[["1"]],"up")))
	}
	return(up)
}

com.lo <- function(x){
	x <- lapply(x,do.call,what=cbind)
	up<-c()
	lo<-c()
	if(any(names(x)=="-1")){
		lo<-rbind(lo,as.matrix(cbind(x[["-1"]],"down")))
	}
	if(any(names(x)=="1")){
		up<-rbind(up,as.matrix(cbind(x[["1"]],"up")))
	}
	return(lo)
}

genes1_ChemVsChem.up<- lapply(comparison1_ChemVsChem_Genes, com.up)
names(genes1_ChemVsChem.up)<-colnames(contrast_matrix)
genes5_ChemVsChem.up<- lapply(comparison5_ChemVsChem_Genes, com.up)
names(genes5_ChemVsChem.up)<-colnames(contrast_matrix)
genes10_ChemVsChem.up<- lapply(comparison10_ChemVsChem_Genes, com.up)
names(genes10_ChemVsChem.up)<-colnames(contrast_matrix)

genes1_ChemVsChem.lo<- lapply(comparison1_ChemVsChem_Genes, com.lo)
names(genes1_ChemVsChem.lo)<-colnames(contrast_matrix)
genes5_ChemVsChem.lo<- lapply(comparison5_ChemVsChem_Genes, com.lo)
names(genes5_ChemVsChem.lo)<-colnames(contrast_matrix)
genes10_ChemVsChem.lo<- lapply(comparison10_ChemVsChem_Genes, com.lo)
names(genes10_ChemVsChem.lo)<-colnames(contrast_matrix)



__________________________________________________________________________________________________________________________

# Chemical classes, dose and KoW analysis

#here we create the factors we are going to use into the anova model and then we do everything as before
colnames(data.ratios)
classes<-as.factor(c(rep("Reactive",times=23),rep("Uncoupler",times=12),rep("Narcotic",times=12),rep("Uncoupler",times=12),rep("Narcotic",times=12),rep("Neurotoxic",times=12),rep("Reactive",times=12),rep("Neurotoxic",times=12),rep("Uncoupler",times=12),rep("Neurotoxic",times=12),rep("Narcotic",times=12)))
KoW<-read.delim("KoW.txt",header=T)
KoW2<-rep(KoW$KoW,each=12)
KoW3<-as.numeric(KoW2[-19])
dose<-rep(c("High","Low"),each=6,times=12)
dose2<-as.factor(dose[-19])
design_aov<-model.matrix(~0+classes+dose2+KoW3)
fit_anova<-lmFit(data.ratios,design_aov)
fit_anova2<-eBayes(fit_anova)

comparison1_Groups<-c()
for (i in 1: length(colnames(design_aov))){
	comparison1_Groups[[i]]<-topTable(fit_anova2, p.value=0.01, number= 50000, genelist=fit_anova2$genes, coef= colnames(design_aov)[[i]], resort.by="t")
	
}

comparison5_Groups<-c()
for (i in 1: length(colnames(design_aov))){
	comparison5_Groups[[i]]<-topTable(fit_anova2, p.value=0.05, number= 50000, genelist=fit_anova2$genes, coef= colnames(design_aov)[[i]], resort.by="t")
	
}

comparison10_Groups<-c()
for (i in 1: length(colnames(design_aov))){
	comparison10_Groups[[i]]<-topTable(fit_anova2, p.value=0.1, number= 50000, genelist=fit_anova2$genes, coef= colnames(design_aov)[[i]], resort.by="t")
	
}

comparison1_Groups.s <- lapply(comparison1_Groups,function(x) try(by(x,sign(x[,"logFC"]),function(y) y)))
comparison1_Groups_Genes<-lapply(comparison1_Groups.s,lapply,function(x){
	if(length(grep("^Error in",x[1])) == 0){
		ann <- annotation.d[match(x$ID,annotation.d$ProbeID),]
	}
	else{
		ann <- NA
	}
	return(list(annotation=ann))
})


comparison5_Groups.s <- lapply(comparison5_Groups,function(x) try(by(x,sign(x[,"logFC"]),function(y) y)))
comparison5_Groups_Genes<-lapply(comparison5_Groups.s,lapply,function(x){
	if(length(grep("^Error in",x[1])) == 0){
		ann <- annotation.d[match(x$ID,annotation.d$ProbeID),]
	}
	else{
		ann <- NA
	}
	return(list(annotation=ann))
})

comparison10_Groups.s <- lapply(comparison10_Groups,function(x) try(by(x,sign(x[,"logFC"]),function(y) y)))
comparison10_Groups_Genes<-lapply(comparison10_Groups.s,lapply,function(x){
	if(length(grep("^Error in",x[1])) == 0){
		ann <- annotation.d[match(x$ID,annotation.d$ProbeID),]
	}
	else{
		ann <- NA
	}
	return(list(annotation=ann))
})

com.up <- function(x){
	x <- lapply(x,do.call,what=cbind)
	up<-c()
	lo<-c()
	if(any(names(x)=="-1")){
		lo<-rbind(lo,as.matrix(cbind(x[["-1"]],"down")))
	}
	if(any(names(x)=="1")){
		up<-rbind(up,as.matrix(cbind(x[["1"]],"up")))
	}
	return(up)
}

com.lo <- function(x){
	x <- lapply(x,do.call,what=cbind)
	up<-c()
	lo<-c()
	if(any(names(x)=="-1")){
		lo<-rbind(lo,as.matrix(cbind(x[["-1"]],"down")))
	}
	if(any(names(x)=="1")){
		up<-rbind(up,as.matrix(cbind(x[["1"]],"up")))
	}
	return(lo)
}

genes1_Groups.up<- lapply(comparison1_Groups_Genes, com.up)
names(genes1_Groups.up)<-colnames(design_aov)
genes5_Groups.up<- lapply(comparison5_Groups_Genes, com.up)
names(genes5_Groups.up)<-colnames(design_aov)
genes10_Groups.up<- lapply(comparison10_Groups_Genes, com.up)
names(genes10_Groups.up)<-colnames(design_aov)

genes1_Groups.lo<- lapply(comparison1_Groups_Genes, com.lo)
names(genes1_Groups.lo)<-colnames(design_aov)
genes5_Groups.lo<- lapply(comparison5_Groups_Genes, com.lo)
names(genes5_Groups.lo)<-colnames(design_aov)
genes10_Groups.lo<- lapply(comparison10_Groups_Genes, com.lo)
names(genes10_Groups.lo)<-colnames(design_aov)



