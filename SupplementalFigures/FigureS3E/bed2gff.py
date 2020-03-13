#!/usr/bin/python3

import sys

etype='enhancer'
if(len(sys.argv) > 3):
	etype=sys.argv[3]

a=0

addregion=0
if len(sys.argv) > 2:
	addregion=int(sys.argv[2])

if sys.argv[1] == '-':
	f=sys.stdin
else:
	f=open(sys.argv[1])

for line in f:
	a+=1
	F=line.strip().split()
	## if we have only coordinates given
	if len(F) < 4:
		F.append('item')

	## gff3 is one base and end included
	print("%s\t%s\t%s\t%i\t%i\t%s=%s" %(F[0],F[3]+'_'+str(a),etype,int(F[1])+1-addregion,int(F[2])+addregion,'\t.\t.\t.\tID',F[3]+'_'+str(a)))
