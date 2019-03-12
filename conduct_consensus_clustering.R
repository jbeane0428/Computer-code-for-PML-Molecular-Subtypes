#load shared library path
shared.library.path <- with(R.version, sprintf("/unprotected/projects/cbmhive/R_packages/R-%s.%s", major, minor));
.libPaths(c(shared.library.path, .libPaths()))

library(biomaRt)
ensembl=useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

library(ConsensusClusterPlus)

#Read in BX discovery data
data<-readRDS(file="outputData/bx.discovery.resid.rds")
dim(data)
gene.set<-readRDS(file="outputData/combined_gene_set.rds")
#subset data to genes in conserved modules
data<-data[unlist(gene.set),]

#normalize data for Consensus Clustering
d.cc = sweep(data,1, apply(data,1,median,na.rm=T))
maxKval<-10

title = "outputData/bx_discovery_resid_ConsensusCluster"
results = ConsensusClusterPlus(d.cc,maxK=maxKval,reps=1000,pItem=0.8,pFeature=1,title=title,clusterAlg="pam",innerLinkage="ward.D",finalLinkage="ward.D",distance="pearson",seed=123456,plot="pdf")
icl = calcICL(results,title=title,plot="pdf")
save(results,icl,file="outputData/bx_discovery_resid_ConsensusCluster.rda")

