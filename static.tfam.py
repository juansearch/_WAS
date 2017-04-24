#!/usr/local/bin/python3
import sys
import re
import os

# cat out.qe862.v2.tfam | ./static.tfam.py QE862.pheno.js.sob.20121216.txt  > QE862.tfam

class Indiv:
	def __init__(self,lst):
		self.fid,self.iid,self.pid,self.mid,self.sex,self.pheno,self.age,self.bmi,self.t2d,self.gender,self.hb1ac,selfq123=lst.split('\t')

indivs={}
head=[]
f1=open(sys.argv[1],'r')
for l in f1:
	if 'pheno' in l:
		head=l.split('\t')
	if 'pheno' not in l:
		i=Indiv(l)
		indivs[i.iid]=i

#fid	iid	pid	mid	sex	pheno	age	bmi	t2d	gender	hb1ac	q123
#0	DGMQ-32087	0	0	1	1	71	16.8	1	1	5.5	2

for l in sys.stdin:
	fid,iid,pid,mid,sex,pheno=l.split(' ')
	sex=indivs[iid].sex
	pheno=indivs[iid].pheno
	print(' '.join([fid,iid,pid,mid,sex,pheno]))

#0	DGMQ-32087	0	0	1	1

