#!/usr/local/bin/python3
import sys
import re
import os

class Iid:
	def __init__(self):
		self.fid=''
		self.iid=''
		self.pid=''
		self.mid=''
		self.gender=''
		self.phenotype=''
		self.age=''
		self.bmi=''
		self.t2d=''
		self.gender=''
		self.hba1c=''
		self.q123=''


db={}
f1=open(sys.argv[1],'r')
for l in f1:
	l=l.rstrip()
	a=l.split('\t')
	if len(a)==12 and 'fid' not in l:
		i=Iid()
		i.fid,i.iid,i.pid,i.mid,i.gender,i.phenotype,i.age,i.bmi,i.t2d,i.gender,i.hba1c,i.q123=a
		db[i.iid]=i
f1.close()


print("fid,iid,gender age,bmi,t2d,pc1,pc2,pc3,pc4,pc5,pc6,pc7,pc8,pc9,pc10")

f2=open(sys.argv[2],'r')
for l in f2:
	l=l.rstrip()
	a=l.split('\t')
	if len(a)==12 and 'FID' not in l:
		if a[1] in db.keys():
			i=db[a[1]]
			i.fid,i.iid,i.pc1,i.pc2,i.pc3,i.pc4,i.pc5,i.pc6,i.pc7,i.pc8,i.pc9,i.pc10=a
			print(','.join([i.fid,i.iid,i.gender,i.age,i.bmi,i.t2d,i.pc1,i.pc2,i.pc3,i.pc4,i.pc5,i.pc6,i.pc7,i.pc8,i.pc9,i.pc10]))
f2.close()
	

#1       2       3       4       5               6       7       8       9       10      11     12
#fid	iid	pid	mid	gender	phenotype	age	bmi	t2d	gender	hba1c	q123

#13      14      15      16      17      18      19      20      21      22      23     24


## add pca to the phenotype file
# print("paste "+staticphe+" "+id+".all.indep.eigenvec | awk -F\"\\t\" '{OFS=\"\\t\";print $13,$14,$10,$7,$8,$9,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}'| perl -pe 's/PC/pc/g' > "+id+".all.pheno.pca.txt")

# paste egaz.v3.pheno.txt t2dg13k.all.indep.eigenvec | awk -F"\t" '{OFS="\t";print $13,$14,$10,$7,$8,$9,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24}'| perl -pe 's/PC/pc/g' > t2dg13k.all.pheno.pca.txt


