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

print("1000g.alt,1000g.ac,1000g.an,1000g.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct,aa_pos_pct,distance,errors_warnings_info")
for l in sys.stdin:
	l=l.rstrip()
	if 'NCHROBS' not in l:
		v=l.split(',')
		if len(v) != 28:
			sys.exit(v)
		k='_'.join([v[5],v[6],v[8],v[9]])
		if k in cpraf:
			print(cpraf[k]+','+l)
		if k not in cpraf:
			print('N,0,0,0,'+l)
#snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct,aa_pos_pct,distance,errors_warnings_info
#c1p65872,G,T,0.04938,1296,1,65872,rs796239852,T,G,./.,3,3,G,upstream_gene_variant,MODIFIER,OR4F5,ENSG00000186092,transcript,ENST00000335137,protein_coding,0.0,c.-3219T>G,,0.0,0.0,0.0,3219

