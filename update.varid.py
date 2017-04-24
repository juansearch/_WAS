#!/usr/local/bin/python3
import sys
import re
import os

for l in sys.stdin:
	l=l.rstrip()
	a=l.split(' ')
	a[1]='c'+a[0]+'p'+a[3]
	print(' '.join(a))
		
#1 rs62635284 0 12783 G A G A G A A A G A G A
