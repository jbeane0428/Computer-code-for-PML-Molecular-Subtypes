#load shared library path
shared.library.path <- with(R.version, sprintf("/unprotected/projects/cbmhive/R_packages/R-%s.%s", major, minor));
.libPaths(c(shared.library.path, .libPaths()))

library(biomaRt)

#Create a list of all modules of co-expressed genes across all datasets
ensembl=useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
m.ensembl=useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")

#Read in all WGCNA results
bx.discovery.wgcna<-readRDS(file="outputData/bx.discovery.wgcna.rds")
br.discovery.wgcna<-readRDS(file="outputData/br.discovery.wgcna.rds")
mouse.wgcna<-readRDS(file="outputData/mouse.wgcna.rds")
tcga.scc.wgcna<-readRDS(file="outputData/tcga.scc.wgcna.rds")

#Creates a list of genes for each module in the WGCNA results
create_gene_list<-function(wgcna.results, label){
        mod.col<-c()
        gene.lab<-c()
	gene.list<-c()
        mod.col<-unique(wgcna.results$modules$Color)
        for(i in 1:length(mod.col)){
                ind<-c()
                ind<-which(mod.col[i]==wgcna.results$modules$Color)
                genes<-c()
                genes<-wgcna.results$modules$Gene[ind]
                #covert IDs to human ENSEMBL IDS
                if(label=="tcga"){
                        ens.q<-getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters=c("hgnc_symbol"),mart=ensembl,values=genes)
                        #Exclude ids that don't map uniquely
                        h.ens.t<-table(ens.q[,1])
                        genes<-names(h.ens.t)[which(h.ens.t==1)]
                }
                if(label=="mouse"){
                        ens.q<-getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"),filters=c("ensembl_gene_id"),mart=m.ensembl,values=genes)
                        #Exclude ids that don't map uniquely
                        h.ens.t<-table(ens.q[,2])
                        genes<-names(h.ens.t)[which(h.ens.t==1)]
                }
                name<-paste(label,mod.col[i],sep="_")
		gene.lab<-c(gene.lab,name)
                gene.list<-c(gene.list,list(genes))
        }
	names(gene.list)<-gene.lab
        return(gene.list)
}

#Compute gene lists for each WGCNA analysis 
bx.list<-create_gene_list(bx.discovery.wgcna,"bx")
br.list<-create_gene_list(br.discovery.wgcna,"br")
mouse.list<-create_gene_list(mouse.wgcna,"mouse")
tcga.list<-create_gene_list(tcga.scc.wgcna,"tcga")

wgcna.gene.list<-c(bx.list,br.list,mouse.list,tcga.list)
saveRDS(wgcna.gene.list,file="outputData/wgcna_gene_list.rds")

#Read in processed data matrices and convert all data to Ensembl Human IDs
bx.discovery.resid<-readRDS(file="outputData/bx.discovery.resid.rds")
br.discovery.resid<-readRDS(file="outputData/br.discovery.resid.rds")
mouse.resid<-readRDS(file="outputData/mouse.resid.rds")
tcga.scc<-readRDS(file="outputData/tcga.scc.resid.rds")

#Convert mouse IDs to human Ensembl IDs
ens<-getBM(attributes=c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"),filters=c("ensembl_gene_id"),mart=m.ensembl,values=rownames(mouse.resid))
h.ens.t<-table(ens[,2])
uni.genes<-names(h.ens.t)[which(h.ens.t==1)]
temp.matrix<-mouse.resid[match(ens[match(uni.genes,ens[,2]),1],rownames(mouse.resid)),]
rownames(temp.matrix)<-uni.genes
mouse.resid<-c()
mouse.resid<-temp.matrix

#Convert TCGA gene symbol to unique Ensembl IDs
ens<-getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters=c("hgnc_symbol"),mart=ensembl,values=rownames(tcga.scc))
h.ens.t<-table(ens[,1])
uni.genes<-names(h.ens.t)[which(h.ens.t==1)]
temp.matrix<-tcga.scc[match(ens[match(uni.genes,ens[,1]),2],rownames(tcga.scc)),]
rownames(temp.matrix)<-uni.genes
tcga.scc<-c()
tcga.scc<-temp.matrix

#Finds correlation between PC1 values for each gene set/module within a dataset
find.correlated.modules<-function(data,gene.set,set.label){
	#Z-score normalize data
	data.z<-c()
	data.z<-t(scale(t(data),center=T,scale=T))
	data.pc1.matrix<-c()
	overlap.length<-c()
	g.length<-c()
	gene.lab<-c()
	gene.lab<-names(gene.set)
	inc.names<-c()
	data.pc1.matrix<-c()
	data.cor.matrix<-c()
	#For each genes set, compute PC1 values across dataset and record gene overlap with original gene set
	for(i in 1:length(gene.lab)){
	        print(i)
       		genes<-c()
	        genes<-unlist(gene.set[i])
       		genes.sub<-intersect(genes,rownames(data.z))
        	if(length(genes.sub)!=0){
        		inc.names<-c(inc.names,names(gene.set[i]))
        		g.length<-c(g.length,length(genes))
        		overlap.length<-c(overlap.length, length(genes.sub))
        		data.sub<-data.z[match(genes.sub,rownames(data.z)),]
        		data.pca<-prcomp(data.sub,center=F,scale=F)
        		data.pc1.matrix<-rbind(data.pc1.matrix,data.pca$rotation[,1])
        	}
	}
	perc<-round((overlap.length/g.length),1)
	rownames(data.pc1.matrix)<-paste(inc.names,overlap.length,g.length,perc,sep="_")

	data.pc1.matrix.orig<-data.pc1.matrix
	data.pc1.matrix<-data.pc1.matrix[which(apply(data.pc1.matrix,1,IQR)!=0),]
	data.cor.matrix<-abs(cor(t(data.pc1.matrix)))
	return(data.cor.matrix)
}

bx.cor.matrix<-find.correlated.modules(bx.discovery.resid, wgcna.gene.list, "bx")
saveRDS(bx.cor.matrix,file="outputData/correlated_modules_bx_discovery.rds")

br.cor.matrix<-find.correlated.modules(br.discovery.resid, wgcna.gene.list, "br")
saveRDS(br.cor.matrix,file="outputData/correlated_modules_br_discovery.rds")

mouse.cor.matrix<-find.correlated.modules(mouse.resid, wgcna.gene.list, "mouse")
saveRDS(mouse.cor.matrix,file="outputData/correlated_modules_mouse.rds")

tcga.cor.matrix<-find.correlated.modules(tcga.scc, wgcna.gene.list, "tcga")
saveRDS(tcga.cor.matrix,file="outputData/correlated_modules_tcga_scc.rds")

#Simplify names of rows and columnms
rownames(br.cor.matrix)<-unlist(lapply(strsplit(rownames(br.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))
colnames(br.cor.matrix)<-unlist(lapply(strsplit(colnames(br.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))
rownames(bx.cor.matrix)<-unlist(lapply(strsplit(rownames(bx.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))
colnames(bx.cor.matrix)<-unlist(lapply(strsplit(colnames(bx.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))
rownames(mouse.cor.matrix)<-unlist(lapply(strsplit(rownames(mouse.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))
colnames(mouse.cor.matrix)<-unlist(lapply(strsplit(colnames(mouse.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))
rownames(tcga.cor.matrix)<-unlist(lapply(strsplit(rownames(tcga.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))
colnames(tcga.cor.matrix)<-unlist(lapply(strsplit(colnames(tcga.cor.matrix),"_"), function(x) paste(x[1],x[2],sep="_")))

#create matrices with the same row/columns - pc1 had to be computed across all 4 for all genesets
com.names<-intersect(intersect(intersect(rownames(br.cor.matrix),rownames(bx.cor.matrix)),rownames(tcga.cor.matrix)),rownames(mouse.cor.matrix))
br.cor.matrix2<-br.cor.matrix[match(com.names,rownames(br.cor.matrix)),match(com.names,colnames(br.cor.matrix))]
bx.cor.matrix2<-bx.cor.matrix[match(com.names,rownames(bx.cor.matrix)),match(com.names,colnames(bx.cor.matrix))]
tcga.cor.matrix2<-tcga.cor.matrix[match(com.names,rownames(tcga.cor.matrix)),match(com.names,colnames(tcga.cor.matrix))]
mouse.cor.matrix2<-mouse.cor.matrix[match(com.names,rownames(mouse.cor.matrix)),match(com.names,colnames(mouse.cor.matrix))]

#find indices of each set of modules
bx.ind<-grep("bx_",rownames(bx.cor.matrix2))
br.ind<-grep("br_",rownames(bx.cor.matrix2))
mouse.ind<-grep("mouse_",rownames(bx.cor.matrix2))
tcga.ind<-grep("tcga_",rownames(bx.cor.matrix2))

#Correlation threshold
cor.cut<-0.85

#Create 0/1 matrices based on correlation cutoff
bx.cor.m<-bx.cor.matrix2
bx.cor.m[which(bx.cor.matrix2<cor.cut)]<-0
bx.cor.m[which(bx.cor.matrix2>=cor.cut)]<-1
br.cor.m<-br.cor.matrix2
br.cor.m[which(br.cor.matrix2<cor.cut)]<-0
br.cor.m[which(br.cor.matrix2>=cor.cut)]<-1
m.cor.m<-mouse.cor.matrix2
m.cor.m[which(mouse.cor.matrix2<cor.cut)]<-0
m.cor.m[which(mouse.cor.matrix2>=cor.cut)]<-1
t.cor.m<-tcga.cor.matrix2
t.cor.m[which(tcga.cor.matrix2<cor.cut)]<-0
t.cor.m[which(tcga.cor.matrix2>=cor.cut)]<-1

sum.matrix<-bx.cor.m+br.cor.m+m.cor.m+t.cor.m
#Subset matrix to just biopsy modules (rows) and other datasets (cols)  because we want to select biopsy modules where the genes were correlated with
#a non-biopsy gene set across all datasets
sum.matrix.sub<-sum.matrix[bx.ind,setdiff(1:ncol(sum.matrix),bx.ind)]

#Print out results
for(i in 1:nrow(sum.matrix.sub)){
        ind<-which(sum.matrix.sub[i,]==4)
        mods<-colnames(sum.matrix.sub)[ind]
	if(length(mods)>0){
		print(i)
		print(mods)
	}
}

#For each of the 9 sets of correlated modules, identify a single gene set
#The single gene set will contain only genes in the biopsy module that overlap with at least one other gene set/module
gs.new<-c()
gs.new.names<-c()
gs.sym<-c()
for(i in 1:nrow(sum.matrix.sub)){
        ind<-which(sum.matrix.sub[i,]==4)
        mods<-colnames(sum.matrix.sub)[ind]
        if(length(mods)>0){
        	genes.all<-c()
        	genes.length<-c()
        	for(j in 1:length(mods)){
	                genes<-unlist(wgcna.gene.list[mods[j]])
        	        genes.num<-length(genes)
                	genes.all<-c(genes.all,genes)
                	genes.length<-c(genes.length,genes.num)
        	}
        	t.genes<-table(genes.all)
        	bx.mod<-i
        	n.inter.genes<-names(t.genes)
        	n.inter.genes.bx<-intersect(n.inter.genes,unlist(wgcna.gene.list[bx.mod]))
        	print(i)
        	print(length(unlist(wgcna.gene.list[bx.mod])))
        	print(length(n.inter.genes))
        	print(length(n.inter.genes.bx))
        	if(length(n.inter.genes.bx)>1){
	                gs.new<-c(gs.new,list(n.inter.genes.bx))
        	        gs.ens<-getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),filters=c("ensembl_gene_id"),mart=ensembl,values=n.inter.genes.bx)
                	gs.new.names<-c(gs.new.names,rownames(sum.matrix.sub)[i])
                	gs.sym<-c(gs.sym,list(gs.ens[which(gs.ens[,2]!=""),2]))
        	}
        }
}
names(gs.new)<-gs.new.names
names(gs.sym)<-gs.new.names
length(gs.new.names)
length(unlist(gs.new))

saveRDS(gs.new,file="outputData/combined_gene_set.rds")
saveRDS(gs.sym,file="outputData/combined_gene_set_symbols.rds")

sessionInfo()

