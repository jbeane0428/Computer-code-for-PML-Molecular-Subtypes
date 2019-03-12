#load shared library path
shared.library.path <- with(R.version, sprintf("/unprotected/projects/cbmhive/R_packages/R-%s.%s", major, minor));
.libPaths(c(shared.library.path, .libPaths()))

library(biomaRt)
ensembl=useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
source("wv_alg.R")

zscore<-function(data){
        z.data<-(data-apply(data,1,mean))/apply(data,1,sd)
        return(z.data)
}

perform.wv<-function(train,train.lab,test,zscore){
        if(zscore==1){
                model<-wv.model(zscore(train),as.factor(train.lab),correction=TRUE)
                predict<-predict.wv(model,zscore(test))
        }
        if(zscore==0){
                model<-wv.model(train,as.factor(train.lab),correction=TRUE)
                predict<-predict.wv(model,test)
        }
        return(predict)
}

predict.smokestatus.persample<-function(test.data,label){
	smoke.sig<-read.table(file="inputData/smoking_signature_genes.txt",sep="\t",header=T)
	smoke.data<-read.table(file="inputData/smoking_signature_data.txt",sep="\t",header=T)
	common.genes<-intersect(rownames(test.data),rownames(smoke.sig))
	
	smoke.sig<-smoke.sig[common.genes,]
	smoke.data<-smoke.data[common.genes,]

	smoke.up<-rownames(smoke.sig)[which(smoke.sig$t>0)]
	smoke.dn<-rownames(smoke.sig)[which(smoke.sig$t<0)]

	sample.lab<-rep(0,length(smoke.data))
	sample.lab[grep("Current",colnames(smoke.data))]<-1

	test.data<-test.data[common.genes,]

	sm.predict.z<-perform.wv(smoke.data, as.factor(sample.lab),test.data,1)

	return(sm.predict.z)
}

bx.discovery.resid<-readRDS(file="outputData/bx.discovery.resid.rds")
bx.discovery.sm<-predict.smokestatus.persample(bx.discovery.resid)
br.discovery.resid<-readRDS(file="outputData/br.discovery.resid.rds")
br.discovery.sm<-predict.smokestatus.persample(br.discovery.resid)
bx.validation.resid<-readRDS(file="outputData/bx.validation.resid.rds")
bx.validation.sm<-predict.smokestatus.persample(bx.validation.resid)
br.validation.resid<-readRDS(file="outputData/br.validation.resid.rds")
br.validation.sm<-predict.smokestatus.persample(br.validation.resid)

#Summarize smoking predictions on a per patient level
all.scores<-c(bx.discovery.sm$scores,br.discovery.sm$scores,bx.validation.sm$scores,br.validation.sm$scores)
sample.data<-read.table(file="inputData/GSE109743_SAMPLES.txt",sep="\t",header=T,row.names=1)
sample.data<-sample.data[names(all.scores),]
names(all.scores)<-sample.data$title

pat<-c()
time<-c()
for(i in 1:length(all.scores)){
        parts<-strsplit(names(all.scores)[i],"-")[[1]]
        pat<-c(pat,parts[3])
        time<-c(time,parts[5])
}
combined<-paste(as.character(pat),as.numeric(time),sep="_")
combined.uni<-unique(combined)

class<-rep(0,length(combined.uni))
for(i in 1:length(combined.uni)){
        ind<-which(combined==combined.uni[i])
        m.scores<-mean(all.scores[ind])
        if(m.scores>0){class[i]<-1}
}

class.all<-rep(0,length(combined))
for(i in 1:length(combined.uni)){
        ind<-which(combined==combined.uni[i])
        class.all[ind]<-class[i]
}
names(class.all)<-names(all.scores)

class.all.status<-rep("Current",length(class.all))
class.all.status[which(class.all==0)]<-"Former"
names(class.all.status)<-names(all.scores)
saveRDS(class.all.status,file="outputData/Genomic_smoking_status.rds")
sessionInfo()
