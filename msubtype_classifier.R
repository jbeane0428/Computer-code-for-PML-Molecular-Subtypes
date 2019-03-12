#load shared library path
shared.library.path <- with(R.version, sprintf("/unprotected/projects/cbmhive/R_packages/R-%s.%s", major, minor));
.libPaths(c(shared.library.path, .libPaths()))

library(biomaRt)
ensembl=useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
library(pamr)
library(heatmap3)
library(sva)

predict.using.pam<-function(genes,class,train){

	pred.genes<-genes
	print(length(pred.genes))

	train.new<-train[match(pred.genes,rownames(train)),]
	train.norm = t(scale(t(train.new),center=T,scale=T))

	d.train<-c()
	d.train<-list(x=as.matrix(train.norm),y=class,geneid=rownames(train.norm))

	set.seed(123456789)
	train<-pamr.train(d.train)

	rownames(d.train$x)<-rownames(train.norm)

	results<-pamr.cv(train, d.train)
	pam.genes<-pamr.listgenes(train, d.train, threshold=0)
	pam.predict.train<-pamr.predict(train,d.train$x,threshold=0,type="class")
	return(list(train=train,predict=pam.predict.train,results=pam.genes))
}

bx.discovery.resid<-readRDS(file="outputData/bx.discovery.resid.rds")
classifier.genes<-readRDS(file="outputData/molecular.subtype.classifier.genes.rds")
load(file="outputData/bx_discovery_resid_ConsensusCluster.rda")
class<-results[[4]]$clrs[[1]]
names(class)<-names(results[[4]]$consensusClass)
class.label<-as.character(class)
class.label[which(class=="#1F78B4")]<-"Proliferative"
class.label[which(class=="#33A02C")]<-"Inflammatory"
class.label[which(class=="#A6CEE3")]<-"Secretory"
class.label[which(class=="#B2DF8A")]<-"Normal-like"

bx.discovery.pam<-predict.using.pam(classifier.genes,class,bx.discovery.resid)
saveRDS(bx.discovery.pam,file="outputData/bx_discovery_msubtype_classifier_results.rds")
#Accuracy of prediction in training set
sum(diag(table(class,bx.discovery.pam$predict)))/sum(table(class,bx.discovery.pam$predict))


#Predict molecular subtype in Validation Set
bx.validation.resid<-readRDS(file="outputData/bx.validation.resid.rds")
common.genes<-intersect(rownames(bx.discovery.resid),rownames(bx.validation.resid))
combo.data<-cbind(bx.discovery.resid[common.genes,],bx.validation.resid[common.genes,])
grp<-c(rep(0,ncol(bx.discovery.resid)),rep(1,ncol(bx.validation.resid)))
com.data<-ComBat(dat=combo.data,batch=grp)
com.data.z<-scale(com.data,center=T,scale=T)
bx.validation.com<-com.data[,which(grp==1)]
bx.validation.com.z<-com.data.z[classifier.genes,which(grp==1)]
d.test<-list(x=as.matrix(bx.validation.com.z),geneid=classifier.genes)
dim(d.test$x)
bx.validation.pam<-pamr.predict(bx.discovery.pam$train,d.test$x,threshold=0,type="class")
saveRDS(bx.validation.pam,file="outputData/bx_validation_msubtype_classifier_results.rds")


#Predict brushes as Proliferative or Not based on 4 gene modules that are concordant between
#biopsies and brushes
gene.ind<-c(grep("bx_darkgrey",names(classifier.genes)),grep("bx_lightgreen",names(classifier.genes)),grep("bx_lightyellow",names(classifier.genes)),grep("bx_magenta",names(classifier.genes)))
brush.classifier.genes<-classifier.genes[gene.ind]

class.Prolif.or.Not<-as.character(class)
class.Prolif.or.Not[which(as.character(class)!="#1F78B4")]<-0


bx.discovery.pam.Prolif.or.Not<-predict.using.pam(brush.classifier.genes,class.Prolif.or.Not,bx.discovery.resid)
#Accuracy of prediction in training set
sum(diag(table(class.Prolif.or.Not,bx.discovery.pam.Prolif.or.Not$predict)))/sum(table(class.Prolif.or.Not,bx.discovery.pam.Prolif.or.Not$predict))
saveRDS(bx.discovery.pam.Prolif.or.Not,file="outputData/bx_discovery_msubtype_classifier_results_Prolif.or.Not.rds")

#Predict Prolif or Not in BX validation
common.genes<-intersect(rownames(bx.discovery.resid),rownames(bx.validation.resid))
combo.data<-cbind(bx.discovery.resid[common.genes,],bx.validation.resid[common.genes,])
grp<-c(rep(0,ncol(bx.discovery.resid)),rep(1,ncol(bx.validation.resid)))
com.data<-ComBat(dat=combo.data,batch=grp)
com.data.z<-scale(com.data,center=T,scale=T)
bx.validation.com<-com.data[,which(grp==1)]
bx.validation.com.z<-com.data.z[brush.classifier.genes,which(grp==1)]
d.test<-list(x=as.matrix(bx.validation.com.z),geneid=classifier.genes)
dim(d.test$x)
bx.validation.pam.Prolif.or.Not<-pamr.predict(bx.discovery.pam.Prolif.or.Not$train,d.test$x,threshold=0,type="class")
saveRDS(bx.validation.pam.Prolif.or.Not,file="outputData/bx_validation_msubtype_classifier_results_Prolif.or.Not.rds")

#Predict Prolif or Not in BR discovery
br.discovery.resid<-readRDS(file="outputData/br.discovery.resid.rds")
common.genes<-intersect(rownames(bx.discovery.resid),rownames(br.discovery.resid))
combo.data<-cbind(bx.discovery.resid[common.genes,],br.discovery.resid[common.genes,])
grp<-c(rep(0,ncol(bx.discovery.resid)),rep(1,ncol(br.discovery.resid)))
com.data<-ComBat(dat=combo.data,batch=grp)
com.data.z<-scale(com.data,center=T,scale=T)
br.discovery.com<-com.data[,which(grp==1)]
br.discovery.com.z<-com.data.z[brush.classifier.genes,which(grp==1)]
d.test<-list(x=as.matrix(br.discovery.com.z),geneid=classifier.genes)
dim(d.test$x)
br.discovery.pam.Prolif.or.Not<-pamr.predict(bx.discovery.pam.Prolif.or.Not$train,d.test$x,threshold=0,type="class")
saveRDS(br.discovery.pam.Prolif.or.Not,file="outputData/br_discovery_msubtype_classifier_results_Prolif.or.Not.rds")

#Predict Prolif or Not in BR validation
br.validation.resid<-readRDS(file="outputData/br.validation.resid.rds")
common.genes<-intersect(rownames(bx.discovery.resid),rownames(br.validation.resid))
combo.data<-cbind(bx.discovery.resid[common.genes,],br.validation.resid[common.genes,])
grp<-c(rep(0,ncol(bx.discovery.resid)),rep(1,ncol(br.validation.resid)))
com.data<-ComBat(dat=combo.data,batch=grp)
com.data.z<-scale(com.data,center=T,scale=T)
br.validation.com<-com.data[,which(grp==1)]
br.validation.com.z<-com.data.z[brush.classifier.genes,which(grp==1)]
d.test<-list(x=as.matrix(br.validation.com.z),geneid=classifier.genes)
dim(d.test$x)
br.validation.pam.Prolif.or.Not<-pamr.predict(bx.discovery.pam.Prolif.or.Not$train,d.test$x,threshold=0,type="class")
saveRDS(br.validation.pam.Prolif.or.Not,file="outputData/br_validation_msubtype_classifier_results_Prolif.or.Not.rds")



