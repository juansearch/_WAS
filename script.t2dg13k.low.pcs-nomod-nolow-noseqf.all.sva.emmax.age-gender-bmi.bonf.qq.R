o<-read.table(file="ft2dg13k.low.pcs-nomod-nolow-noseqf.all.sva.emmax.age-gender-bmi.ps",header=F,sep="\t");
names(o)<-c('cp','stat','pvalue');
o<-o[with(o,order(pvalue)),];
d<-dim(o)[1];
o$rank<-seq(1,d);
o$bonf<-0.05/d;
o$bonf.sig<-o$pvalue<o$bonf;
fn<-"ft2dg13k.low.pcs-nomod-nolow-noseqf.all.sva.emmax.age-gender-bmi.bonf.OUT.txt"
write.table(o,fn,sep="\t",quote=F,row.names=F,col.names=T);
fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), main="all.sva.emmax.low.pcs-nomod-nolow-noseqf",...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }
 l <- 1; fastqq2(pchisq(qchisq(o$pvalue,1)/l,1)); mt<-0.05/length(o$pvalue); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(o$pvalue,1)/l,1)); dev.off();
