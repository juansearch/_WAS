#!/usr/local/bin/python3
import sys
import re
import os

if len(sys.argv) <3:
	sys.exit("./get.sva.for.skat.py <skat.bonf.OUT.post.TRUE.txt> <sva.emmax.age-gender-bmi.bonf.OUT.post.TRUE.txt>")
tgtg=[]
# import skat gene list
f1=open(sys.argv[1],'r')
for l in f1:
	l=l.rstrip()
	a=l.split('\t')
	if 'gene' not in l:
		tgtg.append(a[0])
f1.close()

	#geneh='gene\tpvalue\tnmarkerall\tnmarkertest\tgwrank\tgwbonf\tgwbonfsig\tgwbh\tgwbhsig\tt2dg\tt2drank\tt2dbonf\tt2dbonfsig\tt2dbh\tt2dbhsig'


#      0   1     2       3       4       5          6     7        8     9        10       1           2      3         4     5     6         7         8        9          20     21     22     23       24      25      26      27     28    29    30    31   32  33  34   35       36   37   38  39   40   41  42     43     44      45          46                47
#snph='cp\tstat\tpvalue\tgwrank\tgwbonf\tgwbonfsig\tgwbh\tgwbhsig\tt2dg\tt2drank\tt2dbonf\tt2dbonfsig\tt2dbh\tt2dbhsig\tmafa\tmafu\tnchrobsa\tnchrobsu\tcaddraw\tcaddphred\trefac\taltac\taltaf\texacalt\texacac\texacan\texacaf\tkgalt\tkgac\tkgan\tkgaf\tsnp\ta1\ta2\tmaf\tnchrobs\tchr\tpos\tid\tref\talt\tgt\tannon\tgenen\tallele\tannotation\tannotationimpact\tgenename\tgeneid\tfeaturetype\tfeatureid\ttranscriptbiotype\trank\thgvsc\thgvsp\tcdnapospct\tcdspospct'
# import sva snp list
f2=open(sys.argv[2],'r')
o1=open(sys.argv[1].replace('.txt','.SNPs.txt'),'w')
for l in f2:
	l=l.rstrip()
	a=l.split('\t')
	if 'gene' in l:
		o1.write(l+'\n')
	if 'gene' not in l:
		if a[41] in tgtg:
			o1.write(l+'\n')
f2.close()
o1.close()

#   0   1          2          3          4        5        6        7    8    9       10       11       12         13     14     15     16       17      18      19      20     21    22    23    24  25 26 27  28      29  30  31 32  33  34 35     36     37     38         39                40        41
# a1h,a2h,unaff.a1a1,unaff.a1a2,unaff.a2a2,aff.a1a1,aff.a1a2,aff.a2a2,mafa,mafu,nchrobsa,nchrobsu,cadd.raw,cadd.phred,ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct
#A,G,0,2,288,0,2,574,0.001736,0.003448,1152,580,3.632951,23.2,1728,4,0.0023094688221709007,A,252,119212,2.114e-03,A,55,5008,0.0109824,c1p865584,A,G,0.002309,1732,1,865584,rs148711625,G,A,NA,5,2,A,missense_variant,MODERATE,SAMD11,ENSG00000187634,transcript,ENST00000342066,protein_coding,0.214,c.122G>A,p.Arg41Gln,0.08,0.06
