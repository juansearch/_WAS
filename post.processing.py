#!/usr/local/bin/python3
import sys
import re
import os
import time

#                     0    1     2   3                4   5   6
#./post.processing.py fout.qe866.ac2.pcs-nomod-noseqf.all.sva.emmax.age-gender-bmi.bonf.OUT.txt out.qe866.frq.snpeff.ensembl.kgp.exac.alt.cadd.cc.hwe.csv

if len(sys.argv)<3:
	sys.exit("./post.processing.py <skat/sva txt> <annotation.csv>")


# annotation csv header
#==> fout.qe866.ac2.pcs-nomod-noseqf.csv <==
#  0   1          2          3          4        5        6        7        8        9       10       11       12         13     14     15                         16       17      18      19        20     21    22    23        24        25 26 27       28      29  30     31          32  33  34 35     36     37     38               39                40        41              42           43              44                 45    46       47         48           49          50
#a1h,a2h,unaff.a1a1,unaff.a1a2,unaff.a2a2,aff.a1a1,aff.a1a2,aff.a2a2,    mafa,    mafu,nchrobsa,nchrobsu,cadd.raw,cadd.phred,ref.ac,alt.ac,                    alt.af,exac.alt,exac.ac,exac.an,  exac.af,kg.alt,kg.ac,kg.an,    kg.af,      snp,a1,a2,     maf,nchrobs,chr,   pos,         id,ref,alt,gt,anno.n,gene.n,allele,      annotation,annotation.impact,gene.name,        gene.id,feature.type,     feature.id,transcript.biotype, rank,  hgvs.c,    hgvs.p,cdna.pos.pct,cds.pos.pct
#A  ,G  ,0         ,2         ,288       ,0       ,       2,     574,0.001736,0.003448,    1152,     580,3.632951,      23.2,  1728,     4,     0.0023094688221709007,       A,    252, 119212,2.114e-03,     A,   55, 5008,0.0109824,c1p865584, A, G,0.002309,   1732,  1,865584,rs148711625,  G,  A,NA,     5,     2,     A,missense_variant,         MODERATE,   SAMD11,ENSG00000187634,  transcript,ENST00000342066,    protein_coding,0.214,c.122G>A,p.Arg41Gln,        0.08,       0.06

#

class Runstat:
	def __init__(self,fn):
		self.a=fn.split('.')
		self.model=self.a[5]
		self.kin=self.a[4]
		self.freq=self.a[2]
		self.func=self.a[3]
		self.test=0
		self.nom=0
		self.bh=0
		self.bonf=0
		self.t2d=0
		self.t2dnom=0
		self.t2dbh=0
		self.t2dbonf=0
		self.head="model\tkin\tfreq\tfunc\ttest\tnom\tbh\tbonf\tt2d\tt2dnom\tt2dbh\tt2dbonf\tp.nom.given.t2d\tp.t2d.given.nom"
		self.pt2dgivennom=0.0
		self.pnomgivent2d=0.0
		

	def ostr(self):
		return('\t'.join([self.model,self.kin,self.freq,self.func,str(self.test),str(self.nom),str(self.bh),str(self.bonf),str(self.t2d),str(self.t2dnom),str(self.t2dbh),str(self.t2dbonf),str(self.pt2dgivennom),str(self.pnomgivent2d)]))

class Snp:
	def __init__(self,lst):
		self.a1h,self.a2h,self.unaffa1a1,self.unaffa1a2,self.unaffa2a2,self.affa1a1,self.affa1a2,self.affa2a2,self.mafa,self.mafu,self.nchrobsa,self.nchrobsu,self.caddraw,self.caddphred,self.refac,self.altac,self.altaf,self.exacalt,self.exacac,self.exacan,self.exacaf,self.kgalt,self.kgac,self.kgan,self.kgaf,self.snp,self.a1,self.a2,self.maf,self.nchrobs,self.chr,self.pos,self.id,self.ref,self.alt,self.gt,self.annon,self.genen,self.allele,self.annotation,self.annotationimpact,self.genename,self.geneid,self.featuretype,self.featureid,self.transcriptbiotype,self.rank,self.hgvsc,self.hgvsp,self.cdnapospct,self.cdspospct=lst
		self.cp='0'
		self.stat='0'
		self.pvalue='0'
		self.gwrank='0'
		self.gwbonf='0'
		self.gwbonfsig='FALSE'
		self.gwbh='0'
		self.gwbhsig='FALSE'
		self.t2dg='FALSE'
		self.t2drank='0'
		self.t2dbonf='0'
		self.t2dbonfsig='FALSE'
		self.t2dbh='0'
		self.t2dbhsig='FALSE'
	def ostr(self):
		return('\t'.join([self.cp,self.stat,self.pvalue,self.gwrank,self.gwbonf,self.gwbonfsig,self.gwbh,self.gwbhsig,self.t2dg,self.t2drank,self.t2dbonf,self.t2dbonfsig,self.t2dbh,self.t2dbhsig,self.a1h,self.a2h,self.unaffa1a1,self.unaffa1a2,self.unaffa2a2,self.affa1a1,self.affa1a2,self.affa2a2,self.mafa,self.mafu,self.nchrobsa,self.nchrobsu,self.caddraw,self.caddphred,self.refac,self.altac,self.altaf,self.exacalt,self.exacac,self.exacan,self.exacaf,self.kgalt,self.kgac,self.kgan,self.kgaf,self.snp,self.a1,self.a2,self.maf,self.nchrobs,self.chr,self.pos,self.id,self.ref,self.alt,self.gt,self.annon,self.genen,self.allele,self.annotation,self.annotationimpact,self.genename,self.geneid,self.featuretype,self.featureid,self.transcriptbiotype,self.rank,self.hgvsc,self.hgvsp,self.cdnapospct,self.cdspospct ]))

	#snph='cp\tstat\tpvalue\tgwrank\tgwbonf\tgwbonfsig\tgwbh\tgwbhsig\tt2dg\tt2drank\tt2dbonf\tt2dbonfsig\tt2dbh\tt2dbhsig\ta1h\ta2h\tunaff.a1a1\tunaff.a1a2\tunaff.a2a2\taff.a1a1\taff.a1a2\taff.a2a2\tmafa\tmafu\tnchrobsa\tnchrobsu\tcaddraw\tcaddphred\trefac\taltac\taltaf\texacalt\texacac\texacan\texacaf\tkgalt\tkgac\tkgan\tkgaf\tsnp\ta1\ta2\tmaf\tnchrobs\tchr\tpos\tid\tref\talt\tgt\tannon\tgenen\tallele\tannotation\tannotationimpact\tgenename\tgeneid\tfeaturetype\tfeatureid\ttranscriptbiotype\trank\thgvsc\thgvsp\tcdnapospct\tcdspospct'

class Gen:
	def __init__(self):
		self.gene='0'
		self.pvalue='0'
		self.nmarkerall='0'
		self.nmarkertest='0'
		self.gwrank='0'
		self.gwbonf='0'
		self.gwbonfsig='FALSE'
		self.gwbh='0'
		self.gwbhsig='FALSE'
		self.t2dg='FALSE'
		self.t2drank='0'
		self.t2dbonf='0'
		self.t2dbonfsig='FALSE'
		self.t2dbh='0'
		self.t2dbhsig='FALSE'
	def ostr(self):
		return('\t'.join([self.gene,self.pvalue,self.nmarkerall,self.nmarkertest,self.gwrank,self.gwbonf,self.gwbonfsig,self.gwbh,self.gwbhsig,self.t2dg,self.t2drank,self.t2dbonf,self.t2dbonfsig,self.t2dbh,self.t2dbhsig]))

		# ==> fout.qe867.low.pcs-nomod-noseqf.unrel.skat.null.age-bmi-gender-10pc.bonf.OUT.txt <==
		# SetID	P.value	N.Marker.All	N.Marker.Test	rank	bonf	bonf.sig
		# FOXP4	6.48104226561586e-05	7	7	1	3.87206690931619e-06	FALSE

t0=time.perf_counter()

# list of snp or gene tests in all genes
allt=[]

# skat significant genes gw
gwsigg=[]

# skat significant genes bonf
gwbonfsigg=[]

# skat significant genes bh
gwbhsigg=[]

# skat nominal significant genes gw
gwsigg=[]

# header for regular file
allh=''


# known t2d gene list
t2dg=[]

# list of snp or gene tests in t2d genes
t2dt=[]

# skat significant t2d genes
t2dsigg=[]

#header for t2d file
t2dh=''

# dictionary of snp annotations
snpanno={}

# Read list of t2d genes
f0=open('t2d.genes.txt','r')
for l in f0:
	l=l.rstrip()
	t2dg.append(l)
f0.close()

t1=time.perf_counter()
print("done reading gene list "+str(t1-t0))


# Read list of snps with annotation
f00=open(sys.argv[2],'r')

# Get all the snp annotations
for l in f00:
	l=l.rstrip()
	a=l.split(',')
	if 'nchrobsa' in l:
		allh=a
		continue
	snpanno[a[25]]=Snp(a)
f00.close()

t2=time.perf_counter()
print("done reading annotation "+str(t2-t1))


h=0
# SVA test results processing
if '.sva.' in sys.argv[1]:
	# ==> fout.qe867.low.pcs-nomod-noseqf.kin.sva.emmax.age-gender-bmi.bonf.OUT.txt <==
	# cp	stat	pvalue	rank	bonf	bonf.sig
	# c10p46967660	0.3727269485	7.311639918e-05	1	1.26971228319663e-06	FALSE

	# add snp annotation

	infn=sys.argv[1]
	i1=open(infn,'r')

	n=0
	
	# read in the SNPs and p values
	for l in i1:
		l=l.rstrip()
		if 'stat' not in l:
			d=l.split('\t')
			# Check for SNP in annotation
			if d[0]  in snpanno.keys():
				s=snpanno[d[0]]
				# Add statistical test results to annotation object
				s.cp=d[0]
				s.stat=d[1]
				s.pvalue=d[2]
				s.gwrank=d[3]
				s.gwbonf=d[4]
				s.gwbonfsig=d[5]
				# If known T2D gene, add to this list
				if s.genename in t2dg:
					s.t2dg='TRUE'
					t2dt.append(s)
				# Add to full SNP list
				allt.append(s)
		n=n+1
		if n%1000 == 0:
			print(n)
	i1.close()

	print(str(len(allt))+'\tgenes\t'+str(len(t2dt)))

	#  0   1          2          3          4        5        6        7        8        9       10       11       12         13     14     15                         16       17      18      19        20     21    22    23        24        25 26 27       28      29  30     31          32  33  34 35     36     37     38               39                40        41              42           43              44                 45    46       47         48           49          50
	#a1h,a2h,unaff.a1a1,unaff.a1a2,unaff.a2a2,aff.a1a1,aff.a1a2,aff.a2a2,    mafa,    mafu,nchrobsa,nchrobsu,cadd.raw,cadd.phred,ref.ac,alt.ac,                    alt.af,exac.alt,exac.ac,exac.an,  exac.af,kg.alt,kg.ac,kg.an,    kg.af,      snp,a1,a2,     maf,nchrobs,chr,   pos,         id,ref,alt,gt,anno.n,gene.n,allele,      annotation,annotation.impact,gene.name,        gene.id,feature.type,     feature.id,transcript.biotype, rank,  hgvs.c,    hgvs.p,cdna.pos.pct,cds.pos.pct
	#A  ,G  ,0         ,2         ,288       ,0       ,       2,     574,0.001736,0.003448,    1152,     580,3.632951,      23.2,  1728,     4,     0.0023094688221709007,       A,    252, 119212,2.114e-03,     A,   55, 5008,0.0109824,c1p865584, A, G,0.002309,   1732,  1,865584,rs148711625,  G,  A,NA,     5,     2,     A,missense_variant,         MODERATE,   SAMD11,ENSG00000187634,  transcript,ENST00000342066,    protein_coding,0.214,c.122G>A,p.Arg41Gln,        0.08,       0.06

	# Prepare to write out combined stat / anno
	snph='cp\tstat\tpvalue\tgwrank\tgwbonf\tgwbonfsig\tgwbh\tgwbhsig\tt2dg\tt2drank\tt2dbonf\tt2dbonfsig\tt2dbh\tt2dbhsig\ta1h\ta2h\tunaff.a1a1\tunaff.a1a2\tunaff.a2a2\taff.a1a1\taff.a1a2\taff.a2a2\tmafa\tmafu\tnchrobsa\tnchrobsu\tcaddraw\tcaddphred\trefac\taltac\taltaf\texacalt\texacac\texacan\texacaf\tkgalt\tkgac\tkgan\tkgaf\tsnp\ta1\ta2\tmaf\tnchrobs\tchr\tpos\tid\tref\talt\tgt\tannon\tgenen\tallele\tannotation\tannotationimpact\tgenename\tgeneid\tfeaturetype\tfeatureid\ttranscriptbiotype\trank\thgvsc\thgvsp\tcdnapospct\tcdspospct'

	# Initialize output files
	of1=infn.replace('OUT.txt','OUT.post.txt')
	o1=open(of1,'w')
	o1.write(snph+'\n')
	of2=infn.replace('OUT.txt','OUT.post.TRUE.txt')
	o2=open(of2,'w')
	o2.write(snph+'\n')

	of2n=infn.replace('OUT.txt','OUT.post.NOM.txt')
	o2n=open(of2n,'w')
	o2n.write(snph+'\n')

	ntests=len(allt)
	gwntests=ntests
	alpha=0.05
	print("calculate gw bh")

	# Create stat object
	rs=Runstat(sys.argv[1])
	#self.model,self.kin,self.freq,self.func,self.test,self.nom,self.bh,self.bonf,self.t2d,self.t2dnom,self.t2dbh,self.t2dbonf]))
	of=open(sys.argv[1]+'.stat.txt','w')

	# Calculate BH MTC
	for s in allt:
		# Calculate BH
		s.gwbh=str((float(s.gwrank)/float(ntests))*float(alpha))
		# BH test
		if float(s.pvalue)<float(s.gwbh):
			s.gwbhsig='TRUE'
			rs.bh=rs.bh+1	
		# Write significant snps
		if s.gwbonfsig=='TRUE':
			rs.bonf=rs.bonf+1	
		if s.gwbhsig=='TRUE' or s.gwbonfsig=='TRUE':
			o2.write(s.ostr()+'\n')
		# Write nominal snps
		if float(s.pvalue)<alpha:
			o2n.write(s.ostr()+'\n')
			rs.nom=rs.nom+1
		# Write to full list
		o1.write(s.ostr()+'\n')
		rs.test=rs.test+1	

	o1.close()
	o2.close()
	print('all sva')
	print(rs.head)
	print(rs.ostr())
	
	#cp stat pvalue gwrank gwbonf gwbonfsig gwbh gwbhsig t2dg t2drank t2dgbonf t2dbonfsig 2dbh t2dbhsig

	of3=infn.replace('OUT.txt','OUT.post.t2d.txt')
	o3=open(of3,'w')
	o3.write(snph+'\n')
	of4=infn.replace('OUT.txt','OUT.post.t2d.TRUE.txt')
	o4=open(of4,'w')
	o4.write(snph+'\n')
	of4n=infn.replace('OUT.txt','OUT.post.t2d.NOM.txt')
	o4n=open(of4n,'w')
	o4n.write(snph+'\n')



	ntests=len(t2dt)
	alpha=0.05
	print("calculate t2d bh and bonf")
	rank=0
	for s in t2dt:
		# Calculate BH GW
		s.gwbh=str((float(s.gwrank)/float(gwntests))*float(alpha))
		# GW BH sig
		if float(s.pvalue)<float(s.gwbh):
			s.gwbhsig='TRUE'
			rs.t2dbh=rs.t2dbh+1
		# Calcuate T2D BH
		rank=rank+1
		s.t2drank=str(rank)
		s.t2dbh=str((float(s.t2drank)/float(ntests))*float(alpha))
		# Calculate T2D Bonf
		s.t2dbonf=str(float(alpha)/float(ntests))
		# Test T2D BH
		if float(s.pvalue)<float(s.t2dbh):
			s.t2dbhsig='TRUE'
			rs.t2dbh=rs.t2dbh+1
		# Test T2D Bonf
		if float(s.pvalue)<float(s.t2dbonf):
			s.t2dbonfsig='TRUE'
			rs.t2dbonf=rs.t2dbonf+1
		# Write significant snps
		if s.t2dbhsig=='TRUE' or s.t2dbonfsig=='TRUE':
			o4.write(s.ostr()+'\n')
		# Write nominal snps
		if float(s.pvalue)<alpha:
			o4n.write(s.ostr()+'\n')
			rs.t2dnom=rs.t2dnom+1
		# Write all snps
		o3.write(s.ostr()+'\n')
		rs.t2d=rs.t2d+1

	o3.close()
	o4.close()

	print(rs.head)
	print(rs.ostr())

## SKAT Test results processing
if 'skat' in sys.argv[1]:
	infn=sys.argv[1]
	i1=open(infn,'r')
	rs=Runstat(sys.argv[1])

	n=0
	
#SetID	P.value	N.Marker.All	N.Marker.Test	rank	bonf	bonf.sig
#SNAPC2	7.59460095727249e-05	11	11	1	2.72197724427024e-06	FALSE

#[self.gene,self.pvalue,self.nmarkerall,self.nmarkertest,self.gwrank,self.gwbonf,self.gwbonfsig,self.gwbh,self.gwbhsig,self.t2dg,self.t2drank,self.t2dbonf,self.t2dbonfsig,self.t2dbh,self.t2dbhsig]))
	
	# Get gene test results into object
	for l in i1:
		l=l.rstrip()
		if 'P.value' not in l:
			d=l.split('\t')
			g=Gen()
			g.gene=d[0]
			g.pvalue=d[1]
			g.nmarkerall=d[2]
			g.nmarkertest=d[3]
			g.gwrank=d[4]
			g.gwbonf=d[5]
			g.gwbonfsig=d[6]
			if g.gene in t2dg:
				g.t2dg='TRUE'
				t2dt.append(g)
			allt.append(g)
		n=n+1
		if n%1000 == 0:
			print(n)
	i1.close()

	print(str(len(allt))+'\tgenes\t'+str(len(t2dt))+'\tt2d')
	
	geneh='gene\tpvalue\tnmarkerall\tnmarkertest\tgwrank\tgwbonf\tgwbonfsig\tgwbh\tgwbhsig\tt2dg\tt2drank\tt2dbonf\tt2dbonfsig\tt2dbh\tt2dbhsig'

	sig=0
	t2dsig=0
	of1=infn.replace('OUT.txt','OUT.post.txt')
	o1=open(of1,'w')
	o1.write(geneh+'\n')
	of2=infn.replace('OUT.txt','OUT.post.TRUE.txt')
	o2=open(of2,'w')
	o2.write(geneh+'\n')
	of2n=infn.replace('OUT.txt','OUT.post.NOM.txt')
	o2n=open(of2n,'w')
	o2n.write(geneh+'\n')


	ntests=len(allt)
	gwntests=ntests
	alpha=0.05
	print("calculate gw bh")
	for g in allt:
		rs.test=rs.test+1
		# calculate BH
		g.gwbh=str((float(g.gwrank)/float(ntests))*float(alpha))
		# BH test
		if float(g.pvalue)<float(g.gwbh):
			g.gwbhsig='TRUE'
			rs.bh=rs.bh+1
		# Write all tests
		o1.write(g.ostr()+'\n')
		# Significance test
		if g.gwbhsig=='TRUE' or g.gwbonfsig=='TRUE':
			# Write significant list
			o2.write(g.ostr()+'\n')
			sig=sig+1
			gwsigg.append(g.gene)
		if g.gwbonfsig=='TRUE':
			rs.bonf=rs.bonf+1
		# Nominal test
		if float(g.pvalue)<0.05:
			o2n.write(g.ostr()+'\n')
			rs.nom=rs.nom+1

	o1.close()
	o2.close()
	print(rs.head+'\n')
	print(rs.ostr()+'\n')
	
	#cp stat pvalue gwrank gwbonf gwbonfsig gwbh gwbhsig t2dg t2drank t2dgbonf t2dbonfsig 2dbh t2dbhsig

	of3=infn.replace('OUT.txt','OUT.post.t2d.txt')
	o3=open(of3,'w')
	o3.write(geneh+'\n')
	of4=infn.replace('OUT.txt','OUT.post.t2d.TRUE.txt')
	o4=open(of4,'w')
	o4.write(geneh+'\n')
	of4n=infn.replace('OUT.txt','OUT.post.t2d.NOM.txt')
	o4n=open(of4n,'w')
	o4n.write(geneh+'\n')



	ntests=len(t2dt)
	alpha=0.05
	print("calculate t2d bh and bonf")
	rank=0
	for g in t2dt:
		rs.t2d=rs.t2d+1
		# Calculate BH
		g.gwbh=str((float(g.gwrank)/float(gwntests))*float(alpha))
		# BH test gw
		if float(g.pvalue)<float(g.gwbh):
			g.gwbhsig='TRUE'
		rank=rank+1
		# T2D Multiple test correction
		g.t2drank=str(rank)
		g.t2dbh=str((float(g.t2drank)/float(ntests))*float(alpha))
		g.t2dbonf=str(float(alpha)/float(ntests))
		# T2D BH test
		if float(g.pvalue)<float(g.t2dbh):
			g.t2dbhsig='TRUE'
			rs.t2dbh=rs.t2dbh+1
		if float(g.pvalue)<float(g.t2dbonf):
			g.t2dbonfsig='TRUE'
			rs.t2dbonf=rs.t2dbonf+1
		if float(g.pvalue)<alpha:
			rs.t2dnom=rs.t2dnom+1
		o3.write(g.ostr()+'\n')
		if g.t2dbhsig=='TRUE' or g.t2dbonfsig=='TRUE':
			o4.write(g.ostr()+'\n')
			t2dsig=t2dsig+1
			t2dsigg.append(g.gene)

	o3.close()
	o4.close()

	print(rs.head+'\n')
	print(rs.ostr()+'\n')

#==> fout.qe866.ac2.pcs-nomod-noseqf.all.skat.emmax.age-gender-bmi.bonf.OUT.txt <==
#SetID	P.value	N.Marker.All	N.Marker.Test	rank	bonf	bonf.sig
#ACADS	7.5264105985795e-05	10	10	1	2.79064575542781e-06	FALSE

#==> fout.qe866.ac2.pcs-nomod-noseqf.all.sva.emmax.age-gender-bmi.bonf.OUT.txt <==
#cp	stat	pvalue	rank	bonf	bonf.sig
#c19p18583688	0.1620971357	3.544327705e-07	1	2.41388466459073e-07	FALSE




