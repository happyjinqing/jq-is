library(survival)
file<-read.table("e:/PHD_work/STRING/Hallmark enrichment/survival-subtype(p-2class).txt",header=F,sep="\t")
time<-file[,3]
status<-file[,2]
subtype<-file[,4]
fit<-survfit(Surv(time,status)~subtype,data=file)
plot(fit, lty = 2:3)
legend(2400,.8, c("subtype=1","subtype=2"),lty = 2:3)

library(survminer)
library(ggplot2)
library(ggpubr)
library(magrittr)
require("survival")
fit<-survfit(Surv(time,status)~subtype,data=file)
ggsurvplot(fit,data=file,
main="Survival curves",submain="Based on Kaplan-Meier estimates",
pval=T,risk.table="abs_pct",
risk.table.y.text.col=T,risk.table.y=F,surv_median.line="hv",
ggtheme = theme_light(),xlab="Time in days")
