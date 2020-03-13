#!/usr/bin/python3

import sys

inl=0;

with open(sys.argv[1]) as f,open(sys.argv[1]+'.tsv','w') as g:
	for line in f:
		inl+=1
		ls=line.strip().split()
		if inl == 1:
			llen=len(ls)

		if len(ls) != llen:
			print("line %s has less than %d elements -> skipping" %(inl,llen))
			print("%s\t%s%s"%(ls[0],ls[1],(llen-2)*'\t-1'),file=g)
		else:
			print(line,end='',file=g)

