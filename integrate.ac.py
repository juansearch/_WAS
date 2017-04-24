#!/usr/local/bin/python3

import sys
import re
import os

cpraf={}

f1=open(sys.argv[1],'r')
for l in f1:
	l=l.rstrip()
	if 'CHROM' not in l:
		a=l.split('\t')
		if len(a)==6:
			chrom,pos,nalleles,nchr,a1bc,a2bc=l.split('\t')
			a1b,a1ac=a1bc.split(':')
			a2b,a2ac=a2bc.split(':')
			k='_'.join([chrom,pos,a1b,a2b])
			cpraf[k]=','.join([a1ac,a2ac,str((float(a2ac)/float(nchr)))])
f1.close()	
#1       721694  A       G       5       5008    0.000998403
#1       777318  C       G,T     6,67    5008    0.00119808,0.0133786

print("ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct,aa.pos.pct")
for l in sys.stdin:
	l=l.rstrip()
	if 'NCHROBS' not in l:
		if 'kg' not in l:
			v=l.split(',')
			l1=','.join(v[:34])
			k='_'.join([v[13],v[14],v[16],v[17]])
			if k in cpraf:
				print(cpraf[k]+','+l1)
			if k not in cpraf:
				print('0,0,0,'+l1)
#       0       1       2       3         4        5        6        7   8  9 10  11      12  13  14 15  16  17 18     19     20     21         22                23        24      25           26         27                 28   29     30     31           32          33         34       35                   36
#exac.alt,exac.ac,exac.an,exac.af,1000g.alt,1000g.ac,1000g.an,1000g.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct,aa_pos_pct,distance,errors_warnings_info")


#CHROM   POS     N_ALLELES       N_CHR   {ALLELE:COUNT}
#1       65872   2       1296    T:1232  G:64

