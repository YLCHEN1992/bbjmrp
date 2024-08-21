bbjmrp=function(bbj1,bbj2,title="",ftags="BBJ",fpv=0.05,fr=1e-4,fn=1000,SSRs="SSR1",
fqurt=0.25,fnum=20,expnum=1e+5,methods="median",cuting=1,maxfrq=1,minfrq=0,gpset=F,
gpformat="svg",p1.width=6,p1.height=8,p2.width=8,p2.height=6,
p3.width=6,p3.height=8,dpi=600){

title=paste(unlist(strsplit(bbj1,"\\."))[1],unlist(strsplit(bbj2,"\\."))[1],sep=" Vs. ")

if (!require("ggplot2", quietly = TRUE)){install.packages("ggplot2")}
if (!require("svglite", quietly = TRUE)){install.packages("svglite")}
if (!require("readr", quietly = TRUE)){install.packages("readr")}
library("ggplot2")
library("svglite")
library("readr")

tgtable=function(bbj){
mx=cbind(bbj[,1:5],t(sapply(bbj[,10],function(x){unlist(strsplit(x,":"))[-5]})))
rownames(mx)=1:nrow(mx)
colnames(mx)=c("CHR","BP","SNP","REF","MAF","ES","SE","PV","AF")
mx[,8]=10^(as.numeric(mx[,8])*-1);mx}

pftable=function(mx,fpv=0.05,fr=1e-5,fn=1000,fqurt=0.25,fnum=10){
mx=mx[mx$AF>minfrq&mx$AF<maxfrq,]
mx=mx[mx$PV<fpv,];mx=mx[order(mx$PV),]
pnum=fn;if(nrow(mx)<pnum){pnum=nrow(mx)}
mx=mx[1:pnum,];BPS=round(as.numeric(paste(mx$CHR,mx$BP,sep=""))*fr,0)
mx=mx[!duplicated(BPS),];
nums=round(nrow(mx)*fqurt,0)
if(nums<fnum){nums=fnum}
mx=mx[1:nums,];snps=mx$SNP
list(data=mx,snps=snps)}

fbbj=pftable(tgtable(read.table(gzfile(bbj1),sep="\t")),
fpv=fpv,fr=fr,fn=fn,fqurt=fqurt,fnum=fnum)

bbjx=tgtable(read.table(gzfile(bbj2),sep="\t"))
bbjx=bbjx[bbjx$SNP%in%fbbj$snps,]
insnps=intersect(bbjx$SNP,fbbj$snps)
fbbj$data=fbbj$data[match(insnps,(fbbj$data)$SNP),]
bbjx=bbjx[match(insnps,bbjx$SNP),]

fdata=cbind(fbbj$data,bbjx)[,c(1:18)]
colnames(fdata)=c(paste0("exp_",na.omit(colnames(bbjx))),paste0("out_",na.omit(colnames(bbjx))))
x1=as.numeric(fdata[,6]);x1se=as.numeric(fdata[,7])
x3=as.numeric(fdata[,15]);x3se=as.numeric(fdata[,16])

x3e=x3*(x1/abs(x1));x1e=abs(x1)
SSR1=apply(data.frame(1:length(x1),x1,x3),1,
function(n){kk=summary(lm(x3[-n[1]]~-1+x1[-n[1]],weights=1/(x3se[-n[1]]^2)))$coefficients[1,1]
abs(kk*n[2]-n[3])/((kk^2+1)^0.5)})
k=summary(lm(x3~-1+x1,weights = 1/(x3se^2)))$coefficients[1,1]
SSR2=apply(data.frame(x1,x3),1,function(x){abs(k*x[1]-x[2] )/((k^2+1)^0.5)})
SSR=get(SSRs)
ifelse(methods=="median",{l=which(SSR<=median(SSR));("MR-PRESSO cut=median")},
ifelse(methods=="mean",{l=which(SSR<=mean(SSR));("MR-PRESSO cut=mean")},
ifelse(class(methods)==numeric,{l=which(SSR<=methods);("MR-PRESSO cut=numeric")},
{l=which(SSR<=(mean(SSR)-min(SSR))*cuting);print("MR-PRESSO cut=1/4")})))
SNPtag=fdata$SNP[l]

IVM.noflip=summary(lm(x3[l]~-1+x1[l],weights = 1/(x3se[l]^2)))
IVM.flip=summary(lm(x3e[l]~-1+x1e[l],weights = 1/(x3se[l]^2)))
MR.Egger.noflip=summary(lm(x3[l]~1+x1[l],weights = 1/(x3se[l]^2)))
MR.Egger.flip=summary(lm(x3e[l]~1+x1e[l],weights = 1/(x3se[l]^2)))
source("ftp://systempackage.cn/R/TwoSampleMR.R")
Weighted.median.R=mr_weighted_median(x3[l],x1[l],x3se[l],x1se[l])
Weighted.mode.R=mr_weighted_mode(x3[l],x1[l],x3se[l],x1se[l])
Weighted.median.R.noflip=c(Weighted.median.R$b,Weighted.median.R$se,Weighted.median.R$pval)
Weighted.mode.R.noflip=c(Weighted.mode.R$b,Weighted.mode.R$se,Weighted.mode.R$pval)
Weighted.median.R=mr_weighted_median(x3e[l],x1e[l],x3se[l],x1se[l])
Weighted.mode.R=mr_weighted_mode(x3e[l],x1e[l],x3se[l],x1se[l])
Weighted.median.R.flip=c(Weighted.median.R$b,Weighted.median.R$se,Weighted.median.R$pval)
Weighted.mode.R.flip=c(Weighted.mode.R$b,Weighted.mode.R$se,Weighted.mode.R$pval)
names(Weighted.median.R.noflip)=c("Beta","SE","Pvalues");names(Weighted.mode.R.noflip)=c("Beta","SE","Pvalues")
names(Weighted.median.R.flip)=c("Beta","SE","Pvalues");names(Weighted.mode.R.flip)=c("Beta","SE","Pvalues")
N=expnum;K=length(l);R=summary(lm(x3[l]~-1+x1[l],weights = 1/(x3se[l]^2)))$r.squared
F=((N-K-1)/K )*(R^2/(1-R^2))

mods=data.frame(rbind(IVM.noflip$coefficients[,-3],
IVM.flip$coefficients[,-3],
MR.Egger.noflip$coefficients[2,-3],
MR.Egger.flip$coefficients[2,-3],
Weighted.median.R.noflip,
Weighted.median.R.flip,
Weighted.mode.R.noflip,
Weighted.mode.R.flip))

rownames(mods)=c("IVM.noflips","IVM.flips","MR.Egger.noflips","MR.Egger.flips",
"Weighted.median.noflips","Weighted.median.flips","Weighted.mode.noflips","Weighted.mode.flips")
colnames(mods)=c("Beta","SE","Pvalues")

savedata=fdata[l,]
results=as.data.frame(mods)
inforesults=c(SNPtag,K,"F",F,
"MR.Egger.Intersect.noflip",MR.Egger.noflip$coefficients[1,],
"MR.Egger.Intersect.lip",MR.Egger.flip$coefficients[1,])

write.csv(savedata,paste(ftags,"ori",title,"csv",sep="."))
write.csv(results,paste(ftags,"res",title,"csv",sep="."))
writeLines(inforesults,paste(ftags,"inf",title,"txt",sep="."))}