human_PPI<-read.table("e:/PHD_work/STRING/human_gene_links.txt",header=F,sep="\t")
PPI<-as.matrix(human_PPI)
con<-file("e:/PHD_work/GDC/COAD/somatic mutation/sample-(mutation DEG).txt","r")
line<-readLines(con,n=1)
j=1
while(j<131){
a<-unlist(strsplit(line,"\t"))
ID<-a[1]
geneset<-a[-1]
for(i in 1:length(geneset)){
b1<-which(PPI[,1]==geneset[i])
res<-unique(PPI[b1,])
setwd("e:/PHD_work/STRING/personal PPI network")
name=paste(c(1:130),".txt",sep="")
write.table(res,name[j],sep="\t",col.names=F,row.names=F,quote=F,append=T)
}
j=j+1
line<-readLines(con,n=1)
}