#!/usr/local/bin/python3

import sys
import re
import os

cpraf={}

f1=open(sys.argv[1],'r')
for l in f1:
	l=l.rstrip()
	if 'CHROM' not in l:
		chrom,pos,ref,alt,ac,an,af=l.split('\t')
		if ',' not in alt:
			k='_'.join([chrom,pos,ref,alt])
			cpraf[k]=','.join([alt,ac,an,af])
		if ',' in alt:
			alts=alt.split(',')
			acs=ac.split(',')
			afs=af.split(',')
			for i in range(len(alts)):
				k='_'.join([chrom,pos,ref,alts[i]])
				cpraf[k]=','.join([alts[i],acs[i],an,afs[i]])
f1.close()	
#1       721694  A       G       5       5008    0.000998403
#1       777318  C       G,T     6,67    5008    0.00119808,0.0133786

print("exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct,aa.pos.pct")
for l in sys.stdin:
	l=l.rstrip()
	if 'NCHROBS' not in l:
		if '1000g' not in l:
			v=l.split(',')
			l1=','.join(v[:31])
			k='_'.join([v[9],v[10],v[12],v[13]])
			if k in cpraf:
				print(cpraf[k]+','+l1)
			if k not in cpraf:
				print('N,0,0,0,'+l1)
#        0        1        2        3   4  5  6   7       8   9  10 11  12  13
#1000g.alt,1000g.ac,1000g.an,1000g.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct,aa_pos_pct,distance,errors_warnings_info")

