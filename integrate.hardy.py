#!/usr/local/bin/python3
import sys
import re
import os

class Hardy:
	def __init__(self):
		self.a1a1=''
		self.a1a2=''
		self.a2a2=''
		self.cp=''
		self.group=''
		self.a1=''
		self.a2=''	

# CHR            SNP     TEST   A1   A2                 GENO   O(HET)   E(HET)            P 
#   1      c1p721450      ALL    A    G            1/100/765   0.1155   0.1108       0.3548
#   1      c1p721450      AFF    A    G             1/59/516   0.1024   0.1003            1
#   1      c1p721450    UNAFF    A    G             0/41/249   0.1414   0.1314       0.3774
unaff={}
aff={}

f1=open(sys.argv[1],'r')
for l in f1:
	l=l.rstrip()
	if 'GENO' not in l:
		l=re.sub('^\s+','',l)
		l=re.sub('\s+','\t',l)
		chr,cp,group,a1,a2,geno,ohet,ehet,p=l.split('\t')
		h=Hardy()
		genoa=geno.split('/')
		h.a1=a1
		h.a2=a2
		if len(genoa)==3:
			h.a1a1=genoa[0]
			h.a1a2=genoa[1]
			h.a2a2=genoa[2]
		h.cp=cp
		h.group=group
		if group=='AFF':
			aff[cp]=h
		if group=='UNAFF':
			unaff[cp]=h
f1.close()

#       0          1      2      3      4        5       6       7       8      9    10    11    12  13 14 15  16      17  18  19 20  21  22 23     24
#cadd.raw,cadd.phred,ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct,aa.pos.pct
print("a1h,a2h,unaff.a1a1,unaff.a1a2,unaff.a2a2,aff.a1a1,aff.a1a2,aff.a2a2,mafa,mafu,nchrobsa,nchrobsu,cadd.raw,cadd.phred,ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct")
for l in sys.stdin:
	if 'nchrobs' not in l:
		l=l.rstrip()
		a=l.split(',')
		if a[17] in unaff.keys():
			newstuff=[unaff[a[17]].a1, unaff[a[17]].a2, unaff[a[17]].a1a1, unaff[a[17]].a1a2, unaff[a[17]].a2a2, aff[a[17]].a1a1, aff[a[17]].a1a2, aff[a[17]].a2a2]
			print(','.join(newstuff)+','+l)
		if a[17] not in unaff.keys():
			print('0,0,0,0,0,0,0,0,'+l)
		
		

#0,0,0,N,0,0,0,1000g.alt,1000g.ac,1000g.an,1000g.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct
#1232,64,0.04938271604938271,N,0,0,0,N,0,0,0,c1p65872,G,T,0.04938,1296,1,65872,rs796239852,T,G,./.,3,3,G,upstream_gene_variant,MODIFIER,OR4F5,ENSG00000186092,transcript,ENST00000335137,protein_coding,0.0,c.-3219T>G,,0.0,0.0

#mafa,mafu,nchrobsa,nchrobsu,cadd.raw,cadd.phred,ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct
#0.05295,0.07069,1152,580,-0.073361,1.937,1630,102,0.05889145496535797,A,141,18690,7.544e-03,N , 0, 0, 0,c1p721450, A, G,0.05889,1732, 1,721450,rs79773038, G, A,NA,19, 1, A,downstream_gene_variant,MODIFIER,RP11-206L10.9,ENSG00000237491,transcript,ENST00000434264,lincRNA,0.0,n.*1380G>A,,0.0,0.0
#      0       1    2   3         4     5    6   7                   8 9  10    11        12 13 14 15 16        17 18 19      20   21 22     23         24 25 26 17 28 29 30                      31
