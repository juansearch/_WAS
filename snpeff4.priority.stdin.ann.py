#!/usr/local/bin/python
import sys
import re
import os

impact_convert={'HIGH':4,'MODERATE':3,'LOW':2,'MODIFIER':1,'NONE':0}

#1       865665  rs145442390     G       A       .       .       
# SRC=EXOME
# PL=2
# AF=0.000528541
# ANN=A|splice_donor_variant&splice_region_variant&intron_variant|HIGH|AL645608.1|ENSG00000268179|transcript|ENST00000598827|protein_coding|3/5|c.81+1C>T||||||WARNING_TRANSCRIPT_NO_STOP_CODON
# A|missense_variant|MODERATE|SAMD11|ENSG00000187634|transcript|ENST00000342066|protein_coding|3/14|c.203G>A|p.Arg68Gln|286/2551|203/2046|68/681||
# A|missense_variant|MODERATE|SAMD11|ENSG00000187634|transcript|ENST00000420190|protein_coding|3/7|c.203G>A|p.Arg68Gln|292/626|203/537|68/178||WARNING_TRANSCRIPT_NO_STOP_CODON
# A|missense_variant|MODERATE|SAMD11|ENSG00000187634|transcript|ENST00000437963|protein_coding|3/5|c.203G>A|p.Arg68Gln|263/387|203/327|68/108||WARNING_TRANSCRIPT_NO_STOP_CODON
# A|upstream_gene_variant|MODIFIER|SAMD11|ENSG00000187634|transcript|ENST00000341065|protein_coding||c.-3G>A|||||27|WARNING_TRANSCRIPT_NO_START_CODON
# LOF=(AL645608.1|ENSG00000268179|1|1.00)



class ann:
	def __init__(self,info):
		self.allele,self.annotation,self.annotation_impact,self.gene_name,self.gene_id,self.feature_type,self.feature_id,self.transcript_biotype,self.rank,self.hgvs_c,self.hgvs_p,self.cdna_pos_tot,self.cds_pos_tot,self.aa_pos_tot,self.distance,self.errors_warnings_info=info.split('|')
		self.rank_p=0.0
		self.cdna_pos_p=0.0
		self.cds_pos_p=0.0
		self.aa_pos_p=0.0
		self.impact=impact_convert[self.annotation_impact]
		if '/' in self.rank:
			n,d=self.rank.split('/')
			if float(d) > 0:
				self.rank_p=float(n)/float(d)
		if '/' in self.cdna_pos_tot:
			n,d=self.cdna_pos_tot.split('/')
			if float(d) > 0:
				self.cdna_pos_p=float(n)/float(d)
		if '/' in self.cds_pos_tot:
			n,d=self.cds_pos_tot.split('/')
			if float(d) > 0:
				self.cds_pos_p=float(n)/float(d)
		if '/' in self.aa_pos_tot:
			n,d=self.aa_pos_tot.split('/')
			if float(d) > 0:
				self.aa_pos_p=float(n)/float(d)
	def outstr(self):
		return ','.join([self.allele,self.annotation,self.annotation_impact,self.gene_name,self.gene_id,self.feature_type,self.feature_id,self.transcript_biotype,str(round(self.rank_p,3)),self.hgvs_c,self.hgvs_p,str(round(self.cdna_pos_p,3)),str(round(self.cds_pos_p,3)),str(round(self.aa_pos_p,3)),self.distance])


print "chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct,aa_pos_pct,distance"
for l in sys.stdin:
	if '#' not in l:
		l=l.rstrip()
		vcf=l.split('\t')
		info=''
		for i in vcf[7].split(';'):
			if 'ANN=' in i:
				info=i
		info=info.replace('ANN=','')
		raw_annos=info.split(',')
		annos=[]
		genes=[]
		for i in range(len(raw_annos)):
			annos.append(ann(raw_annos[i]))
			genes.append(annos[-1].gene_name)
		max=annos[0]
		if ',' in vcf[4]:
			vcf[4]=vcf[4].split(',')[0]
		# summary statistics
		print ','.join(vcf[:5])+',NA,'+','.join([str(len(annos)),str(len(set(genes)))])+','+max.outstr()
