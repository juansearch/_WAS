#!/usr/local/bin/python3
import sys
import re
import os

		
#1      c1p65893    A    G     0.001147     0.004566        872        438
freqs={}
f1=open(sys.argv[1],'r')
for l in f1:
	l=l.rstrip()
	if '#' not in l:
		l=re.sub('^\s+','',l)
		l=re.sub('\s+','\t',l)
		chr,snp,a1,a2,mafa,mafu,nchrobsa,nchrobsu=l.split('\t')
		freqs[snp]=','.join([mafa,mafu,nchrobsa,nchrobsu])
f1.close()

#       0          1      2      3      4        5       6       7       8      9    10    11    12  13 14 15  16      17  18  19 20  21  22 23     24
#cadd.raw,cadd.phred,ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct,aa.pos.pct
print("mafa,mafu,nchrobsa,nchrobsu,cadd.raw,cadd.phred,ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct")
for l in sys.stdin:
	if 'nchrobs' not in l:
		l=l.rstrip()
		a=l.split(',')
		if a[13] in freqs.keys():
			print(freqs[a[13]]+','+l)
		if a[13] not in freqs.keys():
			print('0,0,0,0,'+l)
		
		

#0,0,0,N,0,0,0,1000g.alt,1000g.ac,1000g.an,1000g.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct
#1232,64,0.04938271604938271,N,0,0,0,N,0,0,0,c1p65872,G,T,0.04938,1296,1,65872,rs796239852,T,G,./.,3,3,G,upstream_gene_variant,MODIFIER,OR4F5,ENSG00000186092,transcript,ENST00000335137,protein_coding,0.0,c.-3219T>G,,0.0,0.0
