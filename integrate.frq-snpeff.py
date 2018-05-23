#!/usr/local/bin/python3
import sys
import re
import os

# read freq
# t2dg13k.frq.csv <==
# SNP,A1,A2,MAF,NCHROBS
# c1p861292,G,C,0.0007087,25398
frq={}

for l in sys.stdin:
	a=l.rstrip().split(',')
	frq[a[0]]=a


f1=open(sys.argv[1],'r')

# t2dg13k.snpeff.ensembl.csv <==
# chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct,aa_pos_pct,distance

print("SNP,A1,A2,MAF,NCHROBS,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct,aa_pos_pct,distance")
for l in f1:
	a=l.rstrip().split(',')
	cp='c'+a[0]+'p'+a[1]
	if cp in frq.keys():
		print(','.join(frq[cp])+','+','.join(a))

