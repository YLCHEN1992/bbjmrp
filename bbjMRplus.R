bbjMRplus=function(xor){
data=read.csv(xor)
data=data[,-1]
x1=as.numeric(data[,6]);x1se=as.numeric(data[,7])
x3=apply(data[,c(5,14,15)],1,function(x){ifelse(x[1]!=x[2],-1*as.numeric(x[3]),as.numeric(x[3]))});
x3se=as.numeric(data[,16])
x3e=x3*(x1/abs(x1));x1e=abs(x1)

IVMs=as.data.frame(summary(lm(x3e[-1]~-1+x1e[-1],weights = 1/(x3se[-1]^2)))$coefficients)
for(i in 2:length(x1e)){IVMs=rbind(IVMs,as.data.frame(summary(lm(x3e[-i]~-1+x1e[-i],weights = 1/(x3se[-i]^2)))$coefficients))}
metaivm=metafor::rma(yi=IVMs[,1],sei=IVMs[,2],weights = 1/(x3se^2),method='FE')
resuteivmQ=c(metaivm$I2,metaivm$H2,metaivm$QE,metaivm$QEp)

Eggs=as.data.frame(summary(lm(x3e[-1]~1+x1e[-1],weights = 1/(x3se[-1]^2)))$coefficients)
for(i in 2:length(x1e)){Eggs=rbind(Eggs,as.data.frame(summary(lm(x3e[-i]~1+x1e[-i],weights = 1/(x3se[-i]^2)))$coefficients))}
Eggs=Eggs[seq(2,nrow(Eggs),2),]
metaegg=metafor::rma(yi=Eggs[,1],sei=Eggs[,2],weights = 1/(x3se^2),method='EE')
resuteggQ=c(metaegg$I2,metaegg$H2,metaegg$QE,metaegg$QEp)

heterpgeneity=c(paste("IVM(I2,H2,Q,pvalue): ",paste(resuteivmQ,collapse = ",")),paste("MR.Egger(I2,H2,Q,pvalue): ",paste(resuteggQ,collapse = ",")))
writeLines(heterpgeneity,paste0(xor,"mr_heterpgeneity_result.txt"))

IVM.noflip=summary(lm(x3~-1+x1,weights = 1/(x3se^2)))
IVM.flip=summary(lm(x3e~-1+x1e,weights = 1/(x3se^2)))
MR.Egger.noflip=summary(lm(x3~1+x1,weights = 1/(x3se^2)))
MR.Egger.flip=summary(lm(x3e~1+x1e,weights = 1/(x3se^2)))
source("ftp://systempackage.cn/R/TwoSampleMR.R")
Weighted.median.R=mr_weighted_median(x3,x1,x3se,x1se)
Weighted.mode.R=mr_weighted_mode(x3,x1,x3se,x1se)
Weighted.median.R.noflip=c(Weighted.median.R$b,Weighted.median.R$se,Weighted.median.R$pval)
Weighted.mode.R.noflip=c(Weighted.mode.R$b,Weighted.mode.R$se,Weighted.mode.R$pval)
Weighted.median.R=mr_weighted_median(x3e,x1e,x3se,x1se)
Weighted.mode.R=mr_weighted_mode(x3e,x1e,x3se,x1se)
Weighted.median.R.flip=c(Weighted.median.R$b,Weighted.median.R$se,Weighted.median.R$pval)
Weighted.mode.R.flip=c(Weighted.mode.R$b,Weighted.mode.R$se,Weighted.mode.R$pval)
names(Weighted.median.R.noflip)=c("Beta","SE","Pvalues");names(Weighted.mode.R.noflip)=c("Beta","SE","Pvalues")
names(Weighted.median.R.flip)=c("Beta","SE","Pvalues");names(Weighted.mode.R.flip)=c("Beta","SE","Pvalues")
N=441016;K=length(x1);R=summary(lm(x3~-1+x1,weights = 1/(x3se^2)))$r.squared
F=((N-K-1)/K )*(R^2/(1-R^2))
mods=data.frame(rbind(IVM.noflip$coefficients[,-3],
IVM.flip$coefficients[,-3],
MR.Egger.noflip$coefficients[2,-3],
MR.Egger.flip$coefficients[2,-3],
Weighted.median.R.noflip,
Weighted.median.R.flip,
Weighted.mode.R.noflip,
Weighted.mode.R.flip,c(R,K,MR.Egger.flip$coefficients[2,1])))
rownames(mods)=c("IVM.noflips","IVM.flips","MR.Egger.noflips","MR.Egger.flips",
"Weighted.median.noflips","Weighted.median.flips","Weighted.mode.noflips","Weighted.mode.flips","info")
colnames(mods)=c("Beta","SE","Pvalues")

mm=c(MR.Egger.flip$coefficients[2,-3],R,K,MR.Egger.flip$coefficients[1,1],
IVM.flip$coefficients[,4],
Weighted.median.R.flip[3],
Weighted.mode.R.flip[3])
names(mm)=c("Beta","SE","Pvalues","R","K","Int","IVM","Weighted.median","Weighted.mode")
mm}
dd=c()
f=readLines("list")
for(i in 1:length(f)){dd=c(dd,TEST(f[i]))}
mat=matrix(dd,9,length(f))
colnames(mat)=f
rownames(mat)=c("Beta","SE","Pvalues","R","K","Int","IVM","Weighted.median","Weighted.mode")
write.csv(t(mat),"result.csv")