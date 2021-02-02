library(RandomWalkRestartMH)
library(igraph)
i=1
while(i<131){
setwd("e:/PHD_work/STRING/personal PPI network")
name<-paste(c(1:130),".txt",sep="")
PPI_table<-read.table(name[i],sep="\t",header=F)
seedgene<-read.table("colon_cancer_genes(CGC).txt",sep="\t",header=F)
seedgene1<-as.matrix(seedgene)

specific<-PPI_table
PPI_Network <- graph.data.frame(specific,directed=FALSE)
PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)
PPI_MultiplexObject <- create.multiplex(PPI_Network,Layers_Name=c("PPI"))
AdjMatrix_PPI <- compute.adjacency.matrix(PPI_MultiplexObject)
AdjMatrixNorm_PPI <- normalize.multiplex.adjacency(AdjMatrix_PPI)

j=1
setwd("e:/PHD_work/STRING/RWR res0.005")
while(j<62){
n<-which(specific[,1:2]==seedgene1[j])
if(length(n)>0)
{
SeedGene<-seedgene1[j]
RWR_PPI_Results<-Random.Walk.Restart.Multiplex(AdjMatrixNorm_PPI,PPI_MultiplexObject,SeedGene)
result<-RWR_PPI_Results[[1]][which(RWR_PPI_Results[[1]][2]>0.005),]
write.table(result,name[i],quote=F,sep="\t",row.names=F,col.names=F,append=T)
sep=t(c("Gene",SeedGene))
write.table(sep,name[i],quote=F,sep="\t",row.names=F,col.names=F,append=T)
}
j=j+1
}
i=i+1
}



#############
j=1
while(j<131){
setwd("e:/PHD_work/STRING/RWR res0.005")
name<-paste(c(1:130),".txt",sep="")
patient<-read.table(name[j],header=F,sep="\t")
patient1<-as.matrix(patient)
a<-which(patient1[,1]=="Gene")
n1=1
for(i in 1:length(a)){
n2=a[i]-1
if(n2>n1){
geneset<-patient1[n1:n2,1]
res<-c(patient1[a[i],2],geneset)
res1<-t(res)
setwd("e:/PHD_work/STRING/personal gene module0.005")
write.table(res1,name[j],sep="\t",col.names=F,row.names=F,quote=F,append=T)
}
n1=a[i]+1
}
j=j+1
}

##############
j=1
while(j<131){
setwd("e:/PHD_work/STRING/personal gene module0.005")
name<-paste(c(1:130),".txt",sep="")
con<-file(name[j],"r")
line<-readLines(con,n=1)
count<-c()
while(length(line)!=0){
gene<-unlist(strsplit(line,"\t"))
module_size<-length(gene)
res<-c(gene[1],module_size)
res1<-t(res)
setwd("e:/PHD_work/STRING/personal module size0.005")
write.table(res1,name[j],sep="\t",col.names=F,row.names=F,quote=F,append=T)
line<-readLines(con,n=1)
}
j=j+1
}
#############