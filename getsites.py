#!/usr/local/bin/python3
import sys
import re
import os

for l in sys.stdin:
	l=l.rstrip()
	a=l.split('\t')
	if '#CHROM' in l:
		print('\t'.join(a[:8]))
	if '#CHROM' not in l:
		if '#' in l:
			print(l)
	if '#' not in l:
		print('\t'.join(a[:8]))
