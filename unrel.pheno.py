#!/usr/local/bin/python3
import sys
import re
import os

keep=[]
f1=open(sys.argv[1],'r')
#0 DGMQ-31005
for l in f1:
	l=l.rstrip()
	a=l.split(' ')
	keep.append(a[1])
f1.close()

f2=open(sys.argv[2],'r')
#fid	iid	pid	mid	sex	pheno	age	bmi	t2d	gender	hb1ac	q123
#0	DGMQ-32087	0	0	1	1	71	16.8	1	1	5.5	2

for l in f2:
	l=l.rstrip()
	a=l.split('\t')
	if a[1] in keep or 'fid' in l:
		print(l)
f2.close()
	
