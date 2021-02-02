j=1
res1<-c()
res2<-c()
while(j<131){
setwd("e:/PHD_work/STRING/personal PPI network")
name<-paste(c(108:130),".txt",sep="")
patient_PPI<-read.table(name[j],header=F,sep="\t")
PPI<-as.matrix(patient_PPI)
N<-length(unique(c(PPI[,1],PPI[,2])))
setwd("e:/PHD_work/STRING/module integration0.005")
module<-read.table(name[j],header=F,sep="\t")
MG<-as.matrix(module)
n<-length(MG)
con<-file("e:/PHD_work/RWR/hallmark_genes.txt","r")
line<-readLines(con,n=1)
pvalue<-c()
while(length(line)!=0){
a<-unlist(strsplit(line,"\t"))
pathwayID<-a[1]
HG<-a[-1]
M<-length(HG)
k<-length(intersect(HG,MG))
p<-phyper(k-1,M,N-M,n,lower.tail=F)
pvalue<-c(pvalue,p)
line<-readLines(con,n=1)
}
FDR<-p.adjust(pvalue,method='fdr')
res1<-rbind(res1,pvalue)
res2<-rbind(res2,FDR)
j=j+1
}
setwd("e:/PHD_work/STRING")
write.table(res1,"enrichment matrix(p0.005).txt",col.names=F,row.names=F,sep="\t",quote=F,append=T)
write.table(res2,"enrichment matrix(FDR0.005).txt",col.names=F,row.names=F,sep="\t",quote=F,append=T)
