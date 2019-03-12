#load shared library path
shared.library.path <- with(R.version, sprintf("/unprotected/projects/cbmhive/R_packages/R-%s.%s", major, minor));
.libPaths(c(shared.library.path, .libPaths()))

library(heatmap3)
library(edgeR)
library(limma)
library(WGCNA)


#Read in endobronchial biopsy and brush count matrix and sample information from GEO
data<-read.table(file="inputData/GSE109743_exp_count.txt",sep="\t",header=T,row.names=1)
s.data<-read.table(file="inputData/GSE109743_SAMPLES.txt",sep="\t",header=T,row.names=1)

#remove low-quality data
data<-data[,which(is.na(s.data$characteristics..molecular.subtype)==FALSE)]
s.data<-s.data[which(is.na(s.data$characteristics..molecular.subtype)==FALSE),]

#Separate out into biopsy and brush data from discovery and validation cohorts
bx.discovery<-data[,intersect(which(s.data$source.name=="endobronchial biopsy"),which(s.data$characteristics..cohort=="discovery cohort"))]
br.discovery<-data[,intersect(which(s.data$source.name=="normal fluorescing large airway bronchial epithelial cells"),which(s.data$characteristics..cohort=="discovery cohort"))]
bx.validation<-data[,intersect(which(s.data$source.name=="endobronchial biopsy"),which(s.data$characteristics..cohort=="validation cohort"))]
br.validation<-data[,intersect(which(s.data$source.name=="normal fluorescing large airway bronchial epithelial cells"),which(s.data$characteristics..cohort=="validation cohort"))]

compute_residuals<-function(count.matrix, tin, batch){
	data.dge<-c()
	data.dge<-DGEList(counts=count.matrix)
	#Scale Normalization using GLM
	data.dge<-calcNormFactors(data.dge,method="TMM")
	data.cpm<-cpm(data.dge,log=TRUE)
	test1<-apply(data.cpm,1,IQR)
	test2<-apply(data.cpm,1,sum)
	low.genes<-union(which(test1==0),which(test2<=1))
	gene.ind<-setdiff(1:nrow(data.cpm),low.genes)

	data.dge2<-c()
	data.dge2<-data.dge[gene.ind,keep.lib.sizes=F]
	data.dge2<-calcNormFactors(data.dge2,method="TMM")
	data.cpm2<-cpm(data.dge2,log=TRUE)

	design<-model.matrix(~as.numeric(tin)+as.factor(as.character(batch)))
	v<-voom(data.dge2,design=design)
	fit<-lmFit(v,design)
	fit<-eBayes(fit)
	resid<-c()
	resid<-residuals(fit,v)
	return(resid)
}

#Compute residual gene expression values adjusting for tin and batch
bx.discovery.resid<-compute_residuals(bx.discovery,s.data[colnames(bx.discovery),grep("TIN",colnames(s.data))],s.data[colnames(bx.discovery),grep("flowcell.ID",colnames(s.data))])
br.discovery.resid<-compute_residuals(br.discovery,s.data[colnames(br.discovery),grep("TIN",colnames(s.data))],s.data[colnames(br.discovery),grep("flowcell.ID",colnames(s.data))])
bx.validation.resid<-compute_residuals(bx.validation,s.data[colnames(bx.validation),grep("TIN",colnames(s.data))],s.data[colnames(bx.validation),grep("flowcell.ID",colnames(s.data))])
br.validation.resid<-compute_residuals(br.validation,s.data[colnames(br.validation),grep("TIN",colnames(s.data))],s.data[colnames(br.validation),grep("flowcell.ID",colnames(s.data))])

saveRDS(bx.discovery.resid,file="outputData/bx.discovery.resid.rds")
saveRDS(br.discovery.resid,file="outputData/br.discovery.resid.rds")
saveRDS(bx.validation.resid,file="outputData/bx.validation.resid.rds")
saveRDS(br.validation.resid,file="outputData/br.validation.resid.rds")

#Read in mouse count matrix and sample information from GEO
mouse.data<-read.table(file="inputData/GSE111091_rsem_expected_counts.txt",sep="\t",header=T,row.names=1)
mouse.s.data<-read.table(file="inputData/GSE111091_SAMPLES.txt",sep="\t",header=T,row.names=1)
mouse.high.quality.samples<-c("Sample_2B","Sample_40","Sample_41","Sample_42","Sample_45","Sample_46","Sample_47","Sample_48","Sample_49","Sample_5","Sample_7","Sample_70","Sample_71","Sample_72","Sample_73","Sample_74","Sample_75","Sample_76","Sample_78","Sample_79","Sample_7c","Sample_81","Sample_82","Sample_84","Sample_SM1")

#Remove low-quality data
mouse.data<-mouse.data[,match(mouse.high.quality.samples,mouse.s.data$Sample.name)]
mouse.s.data<-mouse.s.data[match(mouse.high.quality.samples,mouse.s.data$Sample.name),]

#Compute residual gene expression values adjusting for tin, strain, and sample type
compute_mouse_residuals<-function(count.matrix, tin, strain, type){
        data.dge<-c()
        data.dge<-DGEList(counts=count.matrix)
        #Scale Normalization using GLM
        data.dge<-calcNormFactors(data.dge,method="TMM")
        data.cpm<-cpm(data.dge,log=TRUE)
        test1<-apply(data.cpm,1,IQR)
        test2<-apply(data.cpm,1,sum)
        low.genes<-union(which(test1==0),which(test2<=1))
        gene.ind<-setdiff(1:nrow(data.cpm),low.genes)

        data.dge2<-c()
        data.dge2<-data.dge[gene.ind,keep.lib.sizes=F]
        data.dge2<-calcNormFactors(data.dge2,method="TMM")
        data.cpm2<-cpm(data.dge2,log=TRUE)

        design<-model.matrix(~as.numeric(tin)+as.factor(as.character(strain))+as.factor(as.character(type)))
        v<-voom(data.dge2,design=design)
        fit<-lmFit(v,design)
        fit<-eBayes(fit)
        resid<-c()
        resid<-residuals(fit,v)
        return(resid)
}

mouse.resid<-compute_mouse_residuals(mouse.data,mouse.s.data$characteristics..median.TIN,mouse.s.data$characteristics..strain,mouse.s.data$source.name)
saveRDS(mouse.resid,file="outputData/mouse.resid.rds")


#Read in TCGA SCC tumor data
tcga.scc.data<-read.table(file="inputData/20150821_LUSC_RNAseqv2_TP_Log2_plus1.txt",sep="\t",header=T)
tcga.scc.annot<-read.table(file="inputData/LUSC_Clinical.txt",sep="\t",header=T)
rownames(tcga.scc.data)<-tcga.scc.data[,1]
tcga.scc.data<-tcga.scc.data[,2:ncol(tcga.scc.data)]
tcga.scc.annot.lab<-gsub("-",".",tcga.scc.annot[,1])

tcga.scc.data.lab<-mapply(function(x) substr(x,1,16),colnames(tcga.scc.data))
tcga.scc.annot.lab2<-mapply(function(x) substr(x,1,16),tcga.scc.annot.lab)

#ordered and subsetted annotation data
tcga.scc.annot.sub<-tcga.scc.annot[match(as.character(tcga.scc.data.lab),as.character(tcga.scc.annot.lab2)),]

#Plate
pl.lab<-c()
pl.lab<-as.character(tcga.scc.annot.sub[,1])
for(i in 1:length(pl.lab)){
        pl.lab[i]<-strtrim(gsub("^TCGA-.{2}-.{4}-.{3}-.{3}-","",pl.lab[i],perl=T),4)
}

tcga.scc.poorquality.samples<-read.table(file="inputData/TCGA_LUSC_PoorQuality_Samples.txt")
tcga.scc.goodsample.ind<-match(setdiff(colnames(tcga.scc.data),tcga.scc.poorquality.samples[,1]),colnames(tcga.scc.data))

#Compute residual gene expression values adjusting for plate
compute_tcga_residuals<-function(cpm.matrix, plate){
	test1<-apply(cpm.matrix,1,IQR)
	test2<-apply(cpm.matrix,1,sum)
	low.genes<-union(which(test1==0),which(test2<=1))
	gene.ind<-setdiff(1:nrow(cpm.matrix),low.genes)

	data.cpm2<-cpm.matrix[gene.ind,]
	design<-model.matrix(~as.factor(plate))
	fit<-lmFit(data.cpm2,design)
	fit<-eBayes(fit)
	resid<-c()
	resid<-residuals(fit,data.cpm2)
        return(resid)
}

tcga.scc.resid<-compute_tcga_residuals(tcga.scc.data[,tcga.scc.goodsample.ind],pl.lab[tcga.scc.goodsample.ind])
saveRDS(tcga.scc.resid,file="outputData/tcga.scc.resid.rds")

#Conduct WGCNA on each dataset
run.wgcna<-function(data,fdr,name){
        source("WGCNA_wrapper.R")
        library(flashClust)
        pdf(file=paste("outputData/",name,"SoftThreshold_plots_fdr",fdr,".pdf",sep=""))
        power <- plotSoftThreshold(data)
        dev.off()
        allowWGCNAThreads()
        # Perform WGCNA
        coexpress <- wgcna(data, softPower=power, saveDissTOM=paste("outputData/",name,"TOM_fdr",fdr,".Rda",sep=""), plotModuleHeatmaps=paste("outputData/",name,"module_heatmaps_fdr",fdr,".pdf",sep=""))
        return(coexpress)
}

bx.discovery.wgcna<-run.wgcna(bx.discovery.resid,0.01,"BX_DISCOVERY")
saveRDS(bx.discovery.wgcna,file="outputData/bx.discovery.wgcna.rds")
br.discovery.wgcna<-run.wgcna(br.discovery.resid,0.01,"BR_DISCOVERY")
saveRDS(br.discovery.wgcna,file="outputData/br.discovery.wgcna.rds")
mouse.wgcna<-run.wgcna(mouse.resid,0.01,"MOUSE")
saveRDS(mouse.wgcna,file="outputData/mouse.wgcna.rds")
tcga.scc.wgcna<-run.wgcna(tcga.scc.resid,0.01,"TCGA_SCC")
saveRDS(tcga.scc.wgcna,file="outputData/tcga.scc.wgcna.rds")

sessionInfo()
