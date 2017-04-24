#!/usr/local/bin/python3
import sys
import re
import os

		
## CADD v1.2 (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013-2015. All rights reserved.
#CHROM	POS	REF	ALT	RawScore	PHRED
#1	65872	T	G	-0.047035	2.688

## CADD v1.2 (c) University of Washington and Hudson-Alpha Institute for Biotechnology 2013-2015. All rights reserved.
#chrom,pos,ref,anc,alt,type,length,istv,isderived,annotype,consequence,consscore,consdetail,gc,cpg,mapability20bp,mapability35bp,scoresegdup,priphcons,mamphcons,verphcons,priphylop,mamphylop,verphylop,gerpn,gerps,gerprs,gerprspval,bstatistic,mutindex,dnahelt,dnamgw,dnaprot,dnaroll,mirsvr-score,mirsvr-e,mirsvr-aln,targetscan,fitcons,chmmtssa,chmmtssaflnk,chmmtxflnk,chmmtx,chmmtxwk,chmmenhg,chmmenh,chmmznfrpts,chmmhet,chmmtssbiv,chmmbivflnk,chmmenhbiv,chmmreprpc,chmmreprpcwk,chmmquies,encexp,ench3k27ac,ench3k4me1,ench3k4me3,encnucleo,encocc,encoccombpval,encocdnasepval,encocfairepval,encocpoliipval,encocctcfpval,encocmycpval,encocdnasesig,encocfairesig,encocpoliisig,encocctcfsig,encocmycsig,segway,toverlapmotifs,motifdist,motifecount,motifename,motifehipos,motifescorechng,tfbs,tfbspeaks,tfbspeaksmax,isknownvariant,esp_af,esp_afr,esp_eur,tg_af,tg_asn,tg_amr,tg_afr,tg_eur,mindisttss,mindisttse,geneid,featureid,ccds,genename,cdnapos,relcdnapos,cdspos,relcdspos,protpos,relprotpos,domain,dst2splice,dst2spltype,exon,intron,oaa,naa,grantham,polyphencat,polyphenval,siftcat,siftval,rawscore,phred

scores={}
f1=open(sys.argv[1],'r')
for l in f1:
	l=l.rstrip()
	if '#' not in l:
		a=l.split('\t')
		if len(a)==6:
			chrom,pos,ref,alt,rawscore,phred=l.split('\t')
			cpra='_'.join([chrom,pos,ref,alt])
			scores[cpra]=rawscore+','+phred
		if len(a) > 6:
			chrom,pos,ref,anc,alt,type,length,istv,isderived,annotype,consequence,consscore,consdetail,gc,cpg,mapabilitytwentybp,mapabilitythirtyfivebp,scoresegdup,priphcons,mamphcons,verphcons,priphylop,mamphylop,verphylop,gerpn,gerps,gerprs,gerprspval,bstatistic,mutindex,dnahelt,dnamgw,dnaprot,dnaroll,mirsvrscore,mirsvre,mirsvraln,targetscan,fitcons,chmmtssa,chmmtssaflnk,chmmtxflnk,chmmtx,chmmtxwk,chmmenhg,chmmenh,chmmznfrpts,chmmhet,chmmtssbiv,chmmbivflnk,chmmenhbiv,chmmreprpc,chmmreprpcwk,chmmquies,encexp,enchthreektwosevenac,enchthreekfourmeone,enchthreekfourmethree,encnucleo,encocc,encoccombpval,encocdnasepval,encocfairepval,encocpoliipval,encocctcfpval,encocmycpval,encocdnasesig,encocfairesig,encocpoliisig,encocctcfsig,encocmycsig,segway,toverlapmotifs,motifdist,motifecount,motifename,motifehipos,motifescorechng,tfbs,tfbspeaks,tfbspeaksmax,isknownvariant,espaf,espafr,espeur,tgaf,tgasn,tgamr,tgafr,tgeur,mindisttss,mindisttse,geneid,featureid,ccds,genename,cdnapos,relcdnapos,cdspos,relcdspos,protpos,relprotpos,domain,dsttwosplice,dsttwospltype,exon,intron,oaa,naa,grantham,polyphencat,polyphenval,siftcat,siftval,rawscore,phred=l.split('\t')
			cpra='_'.join([chrom,pos,ref,alt])
			scores[cpra]=rawscore+','+phred
f1.close()

#     0      1      2        3       4       5       6      7     8     9    10  11 12 13  14      15  16  17 18  19  20 21     22     23     24
#ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct,aa.pos.pct
print("cadd.raw,cadd.phred,ref.ac,alt.ac,alt.af,exac.alt,exac.ac,exac.an,exac.af,kg.alt,kg.ac,kg.an,kg.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno.n,gene.n,allele,annotation,annotation.impact,gene.name,gene.id,feature.type,feature.id,transcript.biotype,rank,hgvs.c,hgvs.p,cdna.pos.pct,cds.pos.pct,aa.pos.pct")
for l in sys.stdin:
	if 'nchrobs' not in l:
		l=l.rstrip()
		a=l.split(',')
		cpra='_'.join([a[16],a[17],a[19],a[20]])
		if cpra in scores.keys():
			print(scores[cpra]+','+l)
		if cpra not in scores.keys():
			print('0,0,'+l)
		
	
		

#0,0,0,N,0,0,0,1000g.alt,1000g.ac,1000g.an,1000g.af,snp,a1,a2,maf,nchrobs,chr,pos,id,ref,alt,gt,anno_n,gene_n,allele,annotation,annotation_impact,gene_name,gene_id,feature_type,feature_id,transcript_biotype,rank,hgvs_c,hgvs_p,cdna_pos_pct,cds_pos_pct
#1232,64,0.04938271604938271,N,0,0,0,N,0,0,0,c1p65872,G,T,0.04938,1296,1,65872,rs796239852,T,G,./.,3,3,G,upstream_gene_variant,MODIFIER,OR4F5,ENSG00000186092,transcript,ENST00000335137,protein_coding,0.0,c.-3219T>G,,0.0,0.0

