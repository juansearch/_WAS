#!/usr/local/bin/python3
import sys
import os
import re

full=0
id0='qe864'
id='qe864'
staticid=id.upper()
phe='pheno.txt'
kgp='./ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz'
exac='./ExAC.r0.3.1.sites.vep.vcf.gz'
ensembl='GRCh37.75'
cadd='./CADD_v1.3/bin/score.sh'
snpeff='java -jar -Xmx40G ./snpEff/snpEff.jar'
vcftools='vcftools'
plink='plink'
bgzip='bgzip'
tabix='tabix'
emmax='./emmax-beta-07Mar2010/emmax'
emmaxkin='./emmax-beta-07Mar2010/emmax-kin'
priority='./snpeff4.priority.stdin.ann.py'
king='./king'
r='./R-3.1.2'
auto='--chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 '
print("date")


class SkatSva:
	def __init__(self):
		self.id=''
		self.filter=''
		self.fn=''

	#f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)
	def skatsva(self,par1,par2,par3,par4,par5,par6,par7,par8,par9,par10,par11,par12):
		# group id
		self.id=par1
		# function filter number
		self.fn=par2
		# frequency and function filter description
		self.filter=par3
		# grep command for functional filter
		self.grep=par4
		# awk command for freq filter
		self.awk=par5
		# file with allele freq
		self.csv=par6
		# file with plink genotypes
		self.plinkbedall=par7
		self.plinkbedun=par8
		# file with phenotypes and PCs
		self.phenoall=par9
		self.phenoun=par10

		# file with kinship matrix
		self.kinfall=par11
		self.kinfun=par12

		self.model=''

		if full==0:
			## SVA KIN

			# filter
			print("head -n 1 "+self.csv+"  > fout."+self.id+"."+self.filter+".csv")
			print("cat "+self.csv+" | "+self.grep+" |  "+self.awk+" > fout."+self.id+"."+self.filter+".csv")
			self.csv2="fout."+self.id+"."+self.filter+".csv"

			## generate setid
			#print("cat "+self.csv2+" | grep -v exac | awk -F\",\" '{print $34,$18}' > fout."+self.id+"."+self.filter+".kin.setid")
			print("cat "+self.csv2+" | grep -v exac | awk -F\",\" '{print $42,$26}' > fout."+self.id+"."+self.filter+".all.setid")

			## filter plink files 
			print("awk '{print $2}' fout."+self.id+"."+self.filter+".all.setid > fout."+self.id+"."+self.filter+".all.setid.keep")
			print(plink+" --bfile "+self.plinkbedall+" --make-bed --out fout."+self.id+"."+self.filter+".all --extract fout."+self.id+"."+self.filter+".all.setid.keep")

			# emmax sva with relatives
			print(plink+" --bfile fout."+self.id+"."+self.filter+".all --recode 12 transpose --out fout."+self.id+"."+self.filter+".all")
			print("cat "+self.phenoall+" | awk -F\"\\"+"t\" '{OFS=\"\\"+"t\";print $1,$2,$6}' | grep -v FID > fout."+self.id+"."+self.filter+".all.sva.emmax.pheno")
			print("cat "+self.phenoall+" | awk -F\"\\"+"t\" '{OFS=\"\\"+"t\";print $1,$2,\"1\",$3,$4,$5,$7,$8}' | grep -v IID > fout."+self.id+"."+self.filter+".all.sva.emmax.cov")
			print(emmax+" -v -d 10  -t fout."+self.id+"."+self.filter+".all -p fout."+self.id+"."+self.filter+".all.sva.emmax.pheno -k "+self.kinfall+" -o fout."+self.id+"."+self.filter+".all.sva.emmax.age-gender-bmi -c fout."+self.id+"."+self.filter+".all.sva.emmax.cov")

			# bonf correct emmax sva with relatives
			self.model='all.sva.emmax.'+self.filter

			# sort, significance, and qq plot of emmax with relatives
			w01=open('script.'+self.id+'.'+self.filter+'.all.sva.emmax.age-gender-bmi.bonf.qq.R','w')
			w01.write("o<-read.table(file=\"fout."+self.id+"."+self.filter+".all.sva.emmax.age-gender-bmi.ps\",header=F,sep=\"\\t\");\n")
			w01.write("names(o)<-c('cp','stat','pvalue');\n")
			w01.write("o<-o[with(o,order(pvalue)),];\n")
			w01.write("d<-dim(o)[1];\n")
			w01.write("o$rank<-seq(1,d);\n")
			w01.write("o$bonf<-0.05/d;\n")
			w01.write("o$bonf.sig<-o$pvalue<o$bonf;\n")
			w01.write("fn<-\"fout."+self.id+"."+self.filter+".all.sva.emmax.age-gender-bmi.bonf.OUT.txt\"\n")
			w01.write("write.table(o,fn,sep=\"\\t\",quote=F,row.names=F,col.names=T);\n")
			w01.write("fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])), main=\""+self.model+"\",...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }\n") 
			w01.write(" l <- 1; fastqq2(pchisq(qchisq(o$pvalue,1)/l,1)); mt<-0.05/length(o$pvalue); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(o$pvalue,1)/l,1)); dev.off();\n")
			w01.close()

			print(r+" CMD BATCH script."+self.id+'.'+self.filter+".all.sva.emmax.age-gender-bmi.bonf.qq.R")
			
			print("./post.processing.py fout."+self.id+'.'+self.filter+".all.sva.emmax.age-gender-bmi.bonf.OUT.txt "+self.csv2)

			svaallpost="fout."+self.id+'.'+self.filter+".all.sva.emmax.age-gender-bmi.bonf.OUT.post.txt"

		## 2 SKAT KIN
		## make r analysis script for skat with kinship matrix


		if full==0:
			self.model='all.skat.emmax.'+self.filter
			w02=open('script.'+self.id+'.'+self.filter+'.all.skat.emmax.age-gender-bmi.bonf.qq.R','w')
			w02.write("library(SKAT);\n")
			w02.write("sessionInfo();\n")
			w02.write("q.bed<-\"./fout."+self.id+"."+self.filter+".all.bed\";\n")
			w02.write("q.bim<-\"./fout."+self.id+"."+self.filter+".all.bim\";\n")
			w02.write("q.fam<-\"./fout."+self.id+"."+self.filter+".all.fam\";\n")
			w02.write("q.setid<-\"./fout."+self.id+"."+self.filter+".all.setid\";\n")
			w02.write("q.ssd<-\"./fout."+self.id+"."+self.filter+".all.ssd\";\n")
			w02.write("q.info<-\"./fout."+self.id+"."+self.filter+".all.info\";\n")
			w02.write("q.cov<-\"./"+self.phenoall+"\";\n")
			w02.write("kf<-\"./"+self.kinfall+"\";\n")
			w02.write("fam_cov<-Read_Plink_FAM_Cov(q.fam,q.cov,Is.binary=FALSE,cov_header=TRUE);\n")
			w02.write("y<-fam_cov$Phenotype;\n")
			w02.write("Generate_SSD_SetID(q.bed,q.bim,q.fam,q.setid,q.ssd,q.info);\n")
			w02.write("ssd.info<-Open_SSD(q.ssd,q.info);\n")
			w02.write("obj<-SKAT_NULL_emmaX(y ~ fam_cov$age + fam_cov$bmi + fam_cov$gender,Kin.File=kf);\n")
			w02.write("out_cov<-SKAT.SSD.All(ssd.info,obj);\n")
			w02.write("t<-out_cov$results;\n")
			w02.write("t<-t[with(t,order(P.value)),];\n")
			w02.write("d<-dim(t)[1];\n")
			w02.write("t$rank<-seq(1,d);\n")
			w02.write("t$bonf<-0.05/d;\n")
			w02.write("t$bonf.sig<-t$P.value<t$bonf;\n")
			w02.write("fn<-\"fout."+self.id+"."+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.txt\";\n")
			w02.write("write.table(t,file=fn,sep=\"\\t\",quote=F,row.names=F,col.names=T);\n")
			w02.write("fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])),  main=\""+self.model+"\",...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }\n") 
			w02.write("l <- 1; fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); mt<-0.05/length(t$P.value); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); dev.off();\n")
			w02.close()
			## run skat emmax with covariates and kinship matrix
			print(r+" CMD BATCH script."+self.id+"."+self.filter+".all.skat.emmax.age-gender-bmi.bonf.qq.R")

			print("./post.processing.py fout."+self.id+'.'+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.txt "+self.csv2)
			skatallpost="fout."+self.id+'.'+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.post.TRUE.txt"
			skatallpost1="fout."+self.id+'.'+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.post.t2d.TRUE.txt"
			skatallpost2="fout."+self.id+'.'+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.post.NOM.txt"
			skatallpost3="fout."+self.id+'.'+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.post.t2d.NOM.txt"

			print("./get.sva.for.skat.py "+skatallpost+" "+svaallpost)
			print("./get.sva.for.skat.py "+skatallpost1+" "+svaallpost)
			print("./get.sva.for.skat.py "+skatallpost2+" "+svaallpost)
			print("./get.sva.for.skat.py "+skatallpost2+" "+svaallpost)

		# UNREL 
		if full==0:

			## generate setid
#			print("cat "+self.csv2+" | grep -v exac | awk -F\",\" '{print $34,$18}' > fout."+self.id+"."+self.filter+".un.setid")
			print("cat "+self.csv2+" | grep -v exac | awk -F\",\" '{print $42,$26}' > fout."+self.id+"."+self.filter+".un.setid")
			print("awk '{print $2}' fout."+self.id+"."+self.filter+".un.setid > fout."+self.id+"."+self.filter+".un.setid.keep")

			## filter plink files 
			print("plink --bfile "+self.plinkbedun+" --make-bed --out fout."+self.id+"."+self.filter+".un --extract fout."+self.id+"."+self.filter+".un.setid.keep")

			# 3 SVA UNREL ONLY WITH KINSHIP MATRIX

			self.model='un.sva.emmax.'+self.filter

			# emmax sva without relatives
			print(plink+" --bfile fout."+self.id+"."+self.filter+".un --recode 12 transpose --out fout."+self.id+"."+self.filter+".un")
			print("cat "+self.phenoun+" | awk -F\"\\"+"t\" '{OFS=\"\\"+"t\";print $1,$2,$6}' | grep -v FID > fout."+self.id+"."+self.filter+".un.sva.emmax.pheno")
			print("cat "+self.phenoun+" | awk -F\"\\"+"t\" '{OFS=\"\\"+"t\";print $1,$2,\"1\",$3,$4,$5,$7,$8}' | grep -v IID > fout."+self.id+"."+self.filter+".un.sva.emmax.cov")
			print(emmax+" -v -d 10  -t fout."+self.id+"."+self.filter+".un -p fout."+self.id+"."+self.filter+".un.sva.emmax.pheno -k "+self.kinfun+" -o fout."+self.id+"."+self.filter+".un.sva.emmax.age-gender-bmi -c fout."+self.id+"."+self.filter+".un.sva.emmax.cov")

			# bonf correct emmax sva with relatives

			# sort, significance, and qq plot of emmax with relatives
			w03=open('script.'+self.id+'.'+self.filter+'.un.sva.emmax.age-gender-bmi.bonf.qq.R','w')
			w03.write("o<-read.table(file=\"fout."+self.id+"."+self.filter+".un.sva.emmax.age-gender-bmi.ps\",header=F,sep=\"\\t\");\n")
			w03.write("names(o)<-c('cp','stat','pvalue');\n")
			w03.write("o<-o[with(o,order(pvalue)),];\n")
			w03.write("d<-dim(o)[1];\n")
			w03.write("o$rank<-seq(1,d);\n")
			w03.write("o$bonf<-0.05/d;\n")
			w03.write("o$bonf.sig<-o$pvalue<o$bonf;\n")
			w03.write("fn<-\"fout."+self.id+"."+self.filter+".un.sva.emmax.age-gender-bmi.bonf.OUT.txt\"\n")
			w03.write("write.table(o,fn,sep=\"\\t\",quote=F,row.names=F,col.names=T);\n")
			w03.write("fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])),  main=\""+self.model+"\",...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }\n") 
			w03.write(" l <- 1; fastqq2(pchisq(qchisq(o$pvalue,1)/l,1)); mt<-0.05/length(o$pvalue); abline(h=-log10(mt), lty=3); png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(o$pvalue,1)/l,1)); dev.off();\n")
			w03.close()

			print(r+" CMD BATCH script."+self.id+'.'+self.filter+".un.sva.emmax.age-gender-bmi.bonf.qq.R")

			print("./post.processing.py fout."+self.id+'.'+self.filter+".un.sva.emmax.age-gender-bmi.bonf.OUT.txt "+self.csv2)
			svaunpost="fout."+self.id+'.'+self.filter+".un.sva.emmax.age-gender-bmi.bonf.OUT.post.txt"

		# 4 SKAT UNREL

		if full==0:
			self.model='un.skat.null.'+self.filter
			# write script
			w04=open('script.'+self.id+'.'+self.filter+'.un.skat.null.age-gender-bmi-10pc.bonf.qq.R','w')
			w04.write("library(SKAT);\n")
			w04.write("sessionInfo();\n")
			w04.write("q.bed<-\"./fout."+self.id+"."+self.filter+".un.bed\";\n")
			w04.write("q.bim<-\"./fout."+self.id+"."+self.filter+".un.bim\";\n")
			w04.write("q.fam<-\"./fout."+self.id+"."+self.filter+".un.fam\";\n")
			w04.write("q.setid<-\"./fout."+self.id+"."+self.filter+".un.setid\";\n")
			w04.write("q.ssd<-\"./fout."+self.id+"."+self.filter+".un.ssd\";\n")
			w04.write("q.info<-\"./fout."+self.id+"."+self.filter+".un.info\";\n")
			w04.write("q.cov<-\"./"+self.phenoun+"\";\n")
			w04.write("fam_cov<-Read_Plink_FAM_Cov(q.fam,q.cov,Is.binary=FALSE,cov_header=TRUE);\n")
			w04.write("y<-fam_cov$Phenotype;\n")
			w04.write("Generate_SSD_SetID(q.bed,q.bim,q.fam,q.setid,q.ssd,q.info);\n")
			w04.write("ssd.info<-Open_SSD(q.ssd,q.info);\n")
			w04.write("obj<-SKAT_Null_Model(y ~ fam_cov$age + fam_cov$bmi + fam_cov$gender + fam_cov$pc1 + fam_cov$pc2 + fam_cov$pc3 + fam_cov$pc4 +fam_cov$pc5 + fam_cov$pc6 + fam_cov$pc7 + fam_cov$pc8 + fam_cov$pc9 + fam_cov$pc10);\n")
			w04.write("out_cov<-SKAT.SSD.All(ssd.info,obj);\n")
			w04.write("t<-out_cov$results;\n")
			w04.write("t<-t[with(t,order(P.value)),];\n")
			w04.write("d<-dim(t)[1];\n")
			w04.write("t$rank<-seq(1,d);\n")
			w04.write("t$bonf<-0.05/d;\n")
			w04.write("t$bonf.sig<-t$P.value<t$bonf;\n")
			w04.write("fn<-\"fout."+self.id+"."+self.filter+".un.skat.null.age-gender-bmi-10pc.bonf.OUT.txt\";\n")
			w04.write("write.table(t,file=fn,sep=\"\\t\",quote=F,row.names=F,col.names=T);\n")
			w04.write("fastqq2 <- function(pvals, ...) { np <- length(pvals); thin.idx <- 1:np; thin.logp.exp <- -log10(thin.idx/(np+1)); thin.logp.obs <- -log10(pvals[order(pvals)[thin.idx]]); plot(thin.logp.exp, thin.logp.obs, xlab=expression(-log[10](p[expected])), ylab=expression(-log[10](p[observed])),  main=\""+self.model+"\",...); abline(0, 1, col='gray', lty=2); thin.idx <- c((0.9)^(5:1), thin.idx); logp.cint.95 <- -log10(qbeta(0.95, thin.idx, np - thin.idx + 1)); logp.cint.05 <- -log10(qbeta(0.05, thin.idx, np - thin.idx + 1)); thin.logp.exp <- -log10(thin.idx/(np+1)); lines(thin.logp.exp, logp.cint.95, lty=2, col='red'); lines(thin.logp.exp, logp.cint.05, lty=2, col='red'); }\n") 
			w04.write("l <- 1;  mt<-0.05/length(t$P.value);  png(paste(fn,'QQ','png',sep='.')); fastqq2(pchisq(qchisq(t$P.value,1)/l,1)); abline(h=-log10(mt), lty=3);dev.off();\n")
			w04.close()

			# run standard skat on unrelateds
			print(r+" CMD BATCH script."+self.id+'.'+self.filter+".un.skat.null.age-gender-bmi-10pc.bonf.qq.R")

			print("./post.processing.py fout."+self.id+'.'+self.filter+".un.skat.null.age-gender-bmi-10pc.bonf.OUT.txt "+self.csv2)
			skatunpost="fout."+self.id+'.'+self.filter+".un.skat.null.age-gender-bmi-10pc.bonf.OUT.post.TRUE.txt"
			skatunpost1="fout."+self.id+'.'+self.filter+".un.skat.null.age-gender-bmi-10pc.bonf.OUT.post.t2d.TRUE.txt"
			skatunpost2="fout."+self.id+'.'+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.post.NOM.txt"
			skatunpost3="fout."+self.id+'.'+self.filter+".all.skat.emmax.age-gender-bmi.bonf.OUT.post.t2d.NOM.txt"
			
			print("./get.sva.for.skat.py "+skatunpost+" "+svaunpost)
			print("./get.sva.for.skat.py "+skatunpost1+" "+svaunpost)
			print("./get.sva.for.skat.py "+skatunpost2+" "+svaunpost)
			print("./get.sva.for.skat.py "+skatunpost3+" "+svaunpost)


####################################################
# pipeline start
if full==0:

	# extract individuals with phenotypes
	print("cat "+staticid+"."+phe+" | grep -v fid | awk '{print $2}' > "+staticid+".keep")
	# exclude singletons, low quality sites, and sites with high missing rate
	print(vcftools+" --gzvcf "+id0+".vcf.gz --keep "+staticid+".keep --mac 2 --max-missing .1 --minQ 30 --recode --out "+id+" "+auto)
	print("mv "+id+".recode.vcf "+id+".vcf")
	print(bgzip+" "+id+".vcf")
	print(tabix+" "+id+".vcf.gz")
	# get site list
	print("zcat "+id+".vcf.gz | ./getsites.py | bgzip -c > "+id+".sites.vcf.gz")
	print(tabix+" "+id+".sites.vcf.gz")


	### VCF to PLINK
	print(plink+" --vcf "+id+".vcf.gz --make-bed --out out."+id+" --const-fid 0")
	# transpose 
	print(plink+" --bfile out."+id+" --recode transpose --out out."+id+".v2")
	# Update TFAM
	print(" cat out."+id+".v2.tfam | ./static.tfam.py "+staticid+"."+phe+"  > "+staticid+".tfam")
	# add variant ids
	print("cat out."+id+".v2.tped | ./update.varid.py > out."+id+".v3.tped")
	print(plink+" --tped out."+id+".v3.tped --tfam "+staticid+".tfam --make-bed --out out."+id+".v4")

	## maximum missing genotype 10% by individual and by site
	print("plink --bfile out."+id+".v4 --make-bed --out out."+id+".v5 --geno 0.1 --mind 0.1 ")


	### POPULATION STRUCTURE FOR OUTLIERS
	# pruning
	print(plink+" --bfile out."+id+" --indep-pairwise 10000 25 0.25 --make-bed --out out."+id+".indep")
	print(plink+" --bfile out."+id+".v5 --indep-pairwise 10000 25 0.25 --make-bed --out out."+id+".v5.indep")

	## get pca
	print(plink+" --bfile out."+id+".indep --pca 10 header tabs --out out."+id+".indep")
	print(plink+" --bfile out."+id+".v5.indep --pca 10 header tabs --out out."+id+".v5.indep")

	# plot PCs
	p01=open('script.out.'+id+'.v5.indep.pca.R','w')
	p01.write('p<-read.table(file="out.'+id+'.indep.eigenvec",header=T);')
	p01.write('names(p)<-c("fid","iid","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10");')
	p01.write('x0<-round(summary(p$pc2)[1],3);')
	p01.write('x1<-round(summary(p$pc2)[6],3);')
	p01.write('xs<-(x1-x0)/10;')
	p01.write('y0<-round(summary(p$pc1)[1],3);')
	p01.write('y1<-round(summary(p$pc1)[6],3);')
	p01.write('ys<-(y1-y0)/10;')
	p01.write('dev.new(height=10,width=10);')
	p01.write('par(mar=c(4,7,1,1));')
	p01.write('pdf("out.'+id+'.indep.pca.pdf");')
	p01.write('plot(p$pc2,p$pc1,frame=F,cex=1.5,pch=21,bg="red",las=1,xlab=NA,xlim=c(x0,x1),xaxp=c(x0,x1,10),ylab=NA,ylim=c(y0,y1),yaxp=c(y0,y1,10));')
	p01.write('mtext(side=1,font=2,line=3,text="PC2");')
	p01.write('mtext(side=2,font=2,line=6,text="PC1");')
	p01.write('axis(side=1,font=2,lwd=2,las=1,at=seq(x0,x1,xs));')
	p01.write('axis(side=2,font=2,lwd=2,las=1,at=seq(y0,y1,ys));')
	p01.write('dev.off();')
	p01.write('p<-read.table(file="out.'+id+'.v5.indep.eigenvec",header=T);')
	p01.write('names(p)<-c("fid","iid","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10");')
	p01.write('x0<-round(summary(p$pc2)[1],3);')
	p01.write('x1<-round(summary(p$pc2)[6],3);')
	p01.write('xs<-(x1-x0)/10;')
	p01.write('y0<-round(summary(p$pc1)[1],3);')
	p01.write('y1<-round(summary(p$pc1)[6],3);')
	p01.write('ys<-(y1-y0)/10;')
	p01.write('dev.new(height=10,width=10);')
	p01.write('par(mar=c(4,7,1,1));')
	p01.write('pdf("out.'+id+'.v5.indep.pca.pdf");')
	p01.write('plot(p$pc2,p$pc1,frame=F,cex=1.5,pch=21,bg="red",las=1,xlab=NA,xlim=c(x0,x1),xaxp=c(x0,x1,10),ylab=NA,ylim=c(y0,y1),yaxp=c(y0,y1,10));')
	p01.write('mtext(side=1,font=2,line=3,text="PC2");')
	p01.write('mtext(side=2,font=2,line=6,text="PC1");')
	p01.write('axis(side=1,font=2,lwd=2,las=1,at=seq(x0,x1,xs));')
	p01.write('axis(side=2,font=2,lwd=2,las=1,at=seq(y0,y1,ys));')
	p01.write('dev.off();')
	p01.close()

	print(r+' CMD BATCH script.out.'+id+'.v5.indep.pca.R')

if full == 0:

	### FREQUENCY
	# get exac and 1000g freq
	print(vcftools+" --gzvcf "+kgp+" --positions-overlap "+id+".vcf.gz --get-INFO AC --get-INFO AN --get-INFO AF --out kgp."+id+"-sites")
	print(vcftools+" --gzvcf "+exac+" --positions-overlap "+id+".vcf.gz --get-INFO AC --get-INFO AN --get-INFO AF --out exac."+id+"-sites")

	# get genotype counts in full cohort
	print(plink+" --bfile out."+id+".v5 --hardy --out out."+id+"")

	# get freq in full cohort
	print(plink+" --bfile out."+id+".v5 --freq --out out."+id+"")

	## get case and control freq
	print(plink+" --bfile out."+id+".v5 --freq case-control --out out."+id+"")

	## get ref/alt ac 
	print(vcftools+" --gzvcf "+id+".vcf.gz --counts --out out."+id+"")

	### ANNOTATION

	if full==0:
		## CADD scores
		print(cadd+" "+id+".vcf.gz  "+id.upper()+".cadd13.tsv.gz")
		print("gunzip  "+id.upper()+".cadd13.tsv.gz")

	# functional annotation 
	print(snpeff+" eff "+ensembl+" "+id+".sites.vcf.gz | bgzip -c > out."+id+".sites.snpeff.ensembl.vcf.gz")

	# get site list
	print("awk '{print $2}' out."+id+".frq | grep -P -v 'SNP' > out."+id+".cp")
	print("cat out."+id+".cp | perl -pe 's/c//g;s/p/\\t/g' > out."+id+".cp.txt")

	# limit annotation to site list
	print(vcftools+" --gzvcf out."+id+".sites.snpeff.ensembl.vcf.gz --positions-overlap out."+id+".cp.txt --recode --recode-INFO-all --out out."+id+".sites.snpeff.ensembl.site-filt")
	print("mv out."+id+".sites.snpeff.ensembl.site-filt.recode.vcf out."+id+".sites.snpeff.ensembl.site-filt.vcf")
	print(bgzip+" out."+id+".sites.snpeff.ensembl.site-filt.vcf")
	print(tabix+" out."+id+".sites.snpeff.ensembl.site-filt.vcf.gz")

	# convert annotation to csv format, one annotation per variant
	print("zcat out."+id+".sites.snpeff.ensembl.site-filt.vcf.gz | "+priority+" > out."+id+".snpeff.ensembl.csv")

	### INTEGRATE
	# frequency and function annotation
	print("cat out."+id+".frq | awk '{OFS=\",\"; print $2,$3,$4,$5,$6}' > out."+id+".frq.csv")
	print("paste out."+id+".frq.csv out."+id+".snpeff.ensembl.csv | perl -pe 's/"+"\\"+"t/,/g' > out."+id+".frq.snpeff.ensembl.csv")
	print("cat out."+id+".frq.snpeff.ensembl.csv | ./integrate.kgp.py kgp."+id+"-sites.INFO | ./integrate.exac.py exac."+id+"-sites.INFO | ./integrate.ac.py out."+id+".frq.count | ./integrate.cadd.py "+id.upper()+".cadd13.tsv | ./integrate.ccfreq.py out."+id+".frq.cc  | ./integrate.hardy.py out."+id+".hwe > out."+id+".frq.snpeff.ensembl.kgp.exac.alt.cadd.cc.hwe.csv")

	## GENERAL ANALYSIS
	## plink files  for general analysis, 
	# plink file of all 
	print(plink+" --bfile out."+id+".v5 --make-bed --out out."+id+".all")

	# transpose for emmax kinship matrix
	print("plink --bfile out."+id+".all --missing-genotype 0 --recode 12 transpose --out out."+id+".all")
	print(emmaxkin+" -v -d 10 out."+id+".all")

	# filter relatives using king
	print("cat out."+id+".all.fam  | awk '{print $2,$2,$3,$4,$5,$6}' > out."+id+".all.fix.fam")
	print(plink+" --bed out."+id+".all.bed --bim out."+id+".all.bim --fam out."+id+".all.fix.fam --make-bed --out out."+id+".all.king-input")
	print(king+" --kinship --ibs -b out."+id+".all.king-input.bed --prefix out."+id+".all.king-out")
	print(king+" -b out."+id+".all.king-input.bed --prefix out."+id+".all.king-out. --unrelated --degree 3")
	print("awk '{print 0,$2}' out."+id+".all.king-out.unrelated.txt > out."+id+".unrel.d3.keep")
	print(plink+" --bfile out."+id+".all --make-bed --keep out."+id+".unrel.d3.keep --out out."+id+".unrel")

	# repeat on unrelated 
	# transpose for emmax kinship matrix
	print(plink+" --bfile out."+id+".unrel  --missing-genotype 0 --recode 12 transpose --out out."+id+".unrel")
	print(emmaxkin+" -v -d 10 out."+id+".unrel")

	### POPULATION STRUCTURE
	# pruning
	print(plink+" --bfile out."+id+".all --indep-pairwise 10000 25 0.25 --make-bed --out out."+id+".all.indep")
	print(plink+" --bfile out."+id+".unrel --indep-pairwise 10000 25 0.25 --make-bed --out out."+id+".unrel.indep")

	## get pca
	print(plink+" --bfile out."+id+".all.indep --pca 10 header tabs --out out."+id+".all.indep")
	print(plink+" --bfile out."+id+".unrel.indep --pca 10 header tabs --out out."+id+".unrel.indep")

	# plot PCs
	p01=open('script.out.'+id+'.indep.pca.R','w')
	p01.write('p<-read.table(file="out.'+id+'.all.indep.eigenvec",header=T);')
	p01.write('names(p)<-c("fid","iid","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10");')
	p01.write('x0<-round(summary(p$pc2)[1],3);')
	p01.write('x1<-round(summary(p$pc2)[6],3);')
	p01.write('xs<-(x1-x0)/10;')
	p01.write('y0<-round(summary(p$pc1)[1],3);')
	p01.write('y1<-round(summary(p$pc1)[6],3);')
	p01.write('ys<-(y1-y0)/10;')
	p01.write('dev.new(height=10,width=10);')
	p01.write('par(mar=c(4,7,1,1));')
	p01.write('pdf("out.'+id+'.all.indep.pca.pdf");')
	p01.write('plot(p$pc2,p$pc1,frame=F,cex=1.5,pch=21,bg="red",las=1,xlab=NA,xlim=c(x0,x1),xaxp=c(x0,x1,10),ylab=NA,ylim=c(y0,y1),yaxp=c(y0,y1,10));')
	p01.write('mtext(side=1,font=2,line=3,text="PC2");')
	p01.write('mtext(side=2,font=2,line=6,text="PC1");')
	p01.write('axis(side=1,font=2,lwd=2,las=1,at=seq(x0,x1,xs));')
	p01.write('axis(side=2,font=2,lwd=2,las=1,at=seq(y0,y1,ys));')
	p01.write('dev.off();')
	p01.write('p<-read.table(file="out.'+id+'.unrel.indep.eigenvec",header=T);')
	p01.write('names(p)<-c("fid","iid","pc1","pc2","pc3","pc4","pc5","pc6","pc7","pc8","pc9","pc10");')
	p01.write('x0<-round(summary(p$pc2)[1],3);')
	p01.write('x1<-round(summary(p$pc2)[6],3);')
	p01.write('xs<-(x1-x0)/10;')
	p01.write('y0<-round(summary(p$pc1)[1],3);')
	p01.write('y1<-round(summary(p$pc1)[6],3);')
	p01.write('ys<-(y1-y0)/10;')
	p01.write('dev.new(height=10,width=10);')
	p01.write('par(mar=c(4,7,1,1));')
	p01.write('pdf("out.'+id+'.unrel.indep.pca.pdf");')
	p01.write('plot(p$pc2,p$pc1,frame=F,cex=1.5,pch=21,bg="red",las=1,xlab=NA,xlim=c(x0,x1),xaxp=c(x0,x1,10),ylab=NA,ylim=c(y0,y1),yaxp=c(y0,y1,10));')
	p01.write('mtext(side=1,font=2,line=3,text="PC2");')
	p01.write('mtext(side=2,font=2,line=6,text="PC1");')
	p01.write('axis(side=1,font=2,lwd=2,las=1,at=seq(x0,x1,xs));')
	p01.write('axis(side=2,font=2,lwd=2,las=1,at=seq(y0,y1,ys));')
	p01.write('dev.off();')
	p01.close()

	print(r+' CMD BATCH script.out.'+id+'.indep.pca.R')

	## add pca to the phenotype file
	print("paste "+staticid+"."+phe+" out."+id+".all.indep.eigenvec | awk -F\"\\t\" '{OFS=\"\\t\";print $13,$14,$10,$7,$8,$9,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}'| perl -pe 's/PC/pc/g' > out."+id+".all.pheno.pca.txt")

	## phenotype file unrelateds only
	print("./unrel.pheno.py out."+id+".unrel.d3.keep "+staticid+"."+phe+" > out."+id+".unrel."+phe)
	print("paste out."+id+".unrel."+phe+" out."+id+".unrel.indep.eigenvec | awk -F\"\\t\" '{OFS=\"\\t\";print $13,$14,$10,$7,$8,$9,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}'| perl -pe 's/PC/pc/g' > out."+id+".unrel.pheno.pca.txt")

	# 1       2       3       4       5       6       7       8       9       10      11      12      13      14      15      16      17      18      19
	# fid	iid	pid	mid	sex	pheno	age	bmi	t2d	gender	hb1ac	q123	FID	IID	PC1	PC2	PC3	PC4	PC5
	# 0	DGMQ-32087	0	0	1	1	71	16.8	1	1	5.5	2	0	DGMQ-32087	-4.73046e-05	-0.027484	0.00803219	1.32471e-05	0.00562024

###############################################################################
### FILTERING

plinkbedall="out."+id+".all"
kinfall="out."+id+".all.BN.kinf"
plinkbedunrel="out."+id+".unrel"
kinfun="out."+id+".unrel.BN.kinf"
phenoall="out."+id+".all.pheno.pca.txt"
phenounrel="out."+id+".unrel.pheno.pca.txt"


csv="out."+id+".frq.snpeff.ensembl.kgp.exac.alt.cadd.cc.hwe.csv"


#1 mafa , 2 mafu , 3 nchrobsa , 4 nchrobsu , 5 cadd.raw , 6 cadd.phred , 7 ref.ac , 8 alt.ac , 9 alt.af , 10 exac.alt , 11 exac.ac , 12 exac.an , 13 exac.af , 14 kg.alt , 15 kg.ac , 16 kg.an , 17 kg.af , 18 snp , 19 a1 , 20 a2 , 21 maf , 22 nchrobs , 23 chr , 24 pos , 25 id , 26 ref , 27 alt , 28 gt , 29 anno.n , 30 gene.n , 31 allele , 32 annotation , 33 annotation.impact , 34 gene.name , 35 gene.id , 36 feature.type , 37 feature.id , 38 transcript.biotype , 39 rank , 40 hgvs.c , 41 hgvs.p , 42 cdna.pos.pct , 43 cds.pos.pct , 
#6 cadd.phred, 7 ref.ac , 8 alt.ac , 9 alt.af 21 maf , 22 nchrobs  

if full==0:
	# FILTER 1 #########################################
	# keep non-singleton variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='1'
	filter="ac2.pcs-nomod-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER' | grep -P -v 'sequence_feature'"
	## Allele frequency filter
	minac='1'
	minmaf='0'
	maxmaf='1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+") {print $0}}'"
	#run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,grep,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 2 #########################################
	# keep rare variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='2'
	filter="rare.pcs-nomod-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER' | grep -P -v 'sequence_feature'"
	## Allele frequency filter
	minac='1'
	minmaf='0'
	maxmaf='.01'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"

	#run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 3 #########################################
	# keep low-freq variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='3'
	filter="low.pcs-nomod-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER' | grep -P -v 'sequence_feature'"

	## Allele frequency filter
	minac='1'
	minmaf='0.01'
	maxmaf='0.1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"

	# run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 4 #########################################
	# keep common variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='4'
	filter="com.pcs-nomod-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER' | grep -P -v 'sequence_feature'"

	## Allele frequency filter
	minac='1'
	minmaf='0.1'
	maxmaf='1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"


	# run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)


	# FILTER 5 #########################################
	# keep potentially deleterious non-singleton variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='5'
	filter="ac2.pcs-nomod-nolow-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER|LOW' | grep -P -v 'sequence_feature'"

	## Allele frequency filter
	minac='1'
	minmaf='0'
	maxmaf='1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"


	#run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 6 #########################################
	# keep potentially deleterious rare variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='6'
	filter="rare.pcs-nomod-nolow-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER|LOW' | grep -P -v 'sequence_feature'"

	## Allele frequency filter
	minac='1'
	minmaf='0'
	maxmaf='.01'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"


	#run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 7 #########################################
	# keep potentially deleterious low-freq variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='7'
	filter="low.pcs-nomod-nolow-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER|LOW' | grep -P -v 'sequence_feature'"

	## Allele frequency filter
	minac='1'
	minmaf='0.01'
	maxmaf='0.1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"


	# run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 8 #########################################
	# keep potentially deleterious common variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='8'
	filter="com.pcs-nomod-nolow-noseqf"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' | grep -P -v 'MODIFIER|LOW' | grep -P -v 'sequence_feature'"

	## Allele frequency filter
	minac='1'
	minmaf='0.1'
	maxmaf='1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"

	# run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 9 #########################################
	# keep all variants in protein coding gene region
	fn='9'
	filter="ac2.pct"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding' "

	## Allele frequency filter
	minac='1'
	minmaf='0'
	maxmaf='1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"


	#run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 10 #########################################
	# keep potentially deleterious rare variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='10'
	filter="rare.pct"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding'"

	## Allele frequency filter
	minac='1'
	minmaf='0'
	maxmaf='.01'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"

	#run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 11 #########################################
	# keep potentially deleterious low-freq variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='11'
	filter="low.pct"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding'"

	## Allele frequency filter
	minac='1'
	minmaf='0.01'
	maxmaf='0.1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"

	# run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	# FILTER 8 #########################################
	# keep potentially deleterious common variants in protein coding sequences in ensembl database assigned by snpeff 
	fn='12'
	filter="com.pct"
	grep="grep -P 'ENSG\d+' | grep -P 'protein_coding'"

	## Allele frequency filter
	minac='1'
	minmaf='0.1'
	maxmaf='1'
	awk=" awk -F\",\" '{if ($15>"+minac+" && $16>"+minac+" && $29>"+minmaf+" && $29<="+maxmaf+") {print $0}}'"

	# run
	f5=SkatSva()
	f5.skatsva(id,fn,filter,grep,awk,csv,plinkbedall,plinkbedunrel,phenoall,phenounrel,kinfall,kinfun)

	print("date")

