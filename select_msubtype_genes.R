#load shared library path
shared.library.path <- with(R.version, sprintf("/unprotected/projects/cbmhive/R_packages/R-%s.%s", major, minor));
.libPaths(c(shared.library.path, .libPaths()))

library(edgeR)
library(limma)
library(biomaRt)
ensembl=useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

#Read in discovery BX data
resid.data<-readRDS(file="outputData/bx.discovery.resid.rds")

#Read in consensus genes
gene.set<-readRDS(file="outputData/combined_gene_set.rds")
#create a vector of genes and their corresponding module colors
mod.colors<-c()
for(i in 1:length(gene.set)){
        mod.colors<-rbind(mod.colors,cbind(unlist(gene.set[i]),rep(names(gene.set)[i],length(unlist(gene.set[i])))))
}
mod.colors[,2]<-gsub("bx_","",mod.colors[,2])

#Load other BR and BX datasets
br.discovery.resid<-readRDS(file="outputData/br.discovery.resid.rds")
br.validation.resid<-readRDS(file="outputData/br.validation.resid.rds")
bx.validation.resid<-readRDS(file="outputData/bx.validation.resid.rds")
common.genes<-intersect(intersect(intersect(rownames(resid.data),rownames(bx.validation.resid)),rownames(br.discovery.resid)),rownames(br.validation.resid))

#Select representative genes for each modules for molecular subtype classifier

#compute eigengenes for each module
d.norm<-t(scale(t(resid.data),center=T,scale=T))
eigen<-c()
for(i in 1:length(gene.set)){
        genes<-unlist(gene.set[i])
        d.pca<-prcomp(d.norm[genes,],center=F,scale=F)
        eigen<-cbind(eigen,d.pca$rotation[,1])
}
cor.pca<-cor(t(resid.data),eigen)^2
cor.max<-apply(cor.pca,1,max)
#Compute the eigengene that each gene is most correlated to
member<-c()
for(i in 1:nrow(d.norm)){
        member<-c(member,which(cor.pca[i,]==cor.max[i]))
}
names(member)<-rownames(d.norm)
#Compute the eigengene assigned to each gene when discovering the consensus modules
orig.member<-c()
for(i in 1:length(gene.set)){
        genes<-unlist(gene.set[i])
        orig.member<-c(orig.member,rep(i,length(genes)))
}

avg.cor<-c()
for(i in 1:length(gene.set)){
        genes<-c()
        mem.genes<-c()
	#Select genes in modules that are in all datasets
        genes<-intersect(unlist(gene.set[i]),common.genes)
        #Find genes that in are in module i by eigengene correlation above
	mem.sub<-member[match(genes,names(member))]
        mem.sub<-mem.sub[which(mem.sub==i)]
        mem.genes<-names(mem.sub)
        print(length(genes))
        print(length(mem.genes))
        d.pca<-prcomp(d.norm[mem.genes,],center=F,scale=F)
        avg.cor<-c(avg.cor,mean(cor(t(resid.data[mem.genes,]),d.pca$rotation[,1])^2))
}
#Include more genes in classifier for each module if there is more noise or
#A lower average correlation to the module eigengene
un.cor<-(1-avg.cor)
final.perc<-(un.cor)/sum(un.cor)
num.genes.for.each.module<-c(round(final.perc*20))


classifier.genes<-c()
classifier.gene.colors<-c()
for(i in 1:length(gene.set)){
        genes<-c()
        mem.genes<-c()
	#Select genes in modules that are in all datasets
        genes<-intersect(unlist(gene.set[i]),common.genes)
	#Find genes that in are in module i by eigengene correlation above
        mem.sub<-member[match(genes,names(member))]
        mem.sub<-mem.sub[which(mem.sub==i)]
        mem.genes<-names(mem.sub)
        print(length(genes))
        print(length(mem.genes))
	#Compute PCA across mem.genes (genes that are in module i by eigengene correlation)
        d.pca<-prcomp(d.norm[mem.genes,],center=F,scale=F)
	#Compute correlation between PC1 and the data
        temp.cor<-c()
        temp.cor<-cor(t(resid.data[mem.genes,]),d.pca$rotation[,1])^2
        names(temp.cor)<-mem.genes
	#Select most correlated genes
        s.temp.cor<-sort(temp.cor,decreasing=T,index.return=T)
        print(length(temp.cor))
        classifier.genes<-c(classifier.genes,names(temp.cor)[s.temp.cor$ix[1:num.genes.for.each.module[i]]])
        classifier.gene.colors<-c(classifier.gene.colors,rep(names(gene.set[i]),num.genes.for.each.module[i]))
}

names(classifier.genes)<-classifier.gene.colors
saveRDS(classifier.genes,file="outputData/molecular.subtype.classifier.genes.rds")
