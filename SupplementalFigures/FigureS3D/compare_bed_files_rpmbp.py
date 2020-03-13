#!/usr/bin/env python3.7

# compare_bed_files_rpmbp.py -  python3 script
# Copyright (C) 2018  Sebastian Mackowiak
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import sys
import os
import subprocess
from subprocess import PIPE,run
import numpy as np
import time
from optparse import OptionParser

import numpy 
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import pandas as pd

def perc(x):
	return int(100*x)

def redof(f,v):
	if not os.path.isfile(f) or v == 1:
		return 1
	else:
		return 0

def fpath(f):
	return os.popen("readlink -f "+f).read().strip()

def intersectBed(k1,k2,p,r=0):
	bd={'wo':'intersect','v':'uniq_'+myd[k1]['ab']}
	cmd='intersectBed -a {0} -b {1} -{2}'.format(myd[k1]['file'],myd[k2]['file'],p)
	if r == 1:
		k1,k2 = k2,k1

	of=myd[k1]['ab']+'_vs_'+myd[k2]['ab']+'.'+bd[p]
	if redof(of,redo):
		with open(of,'w') as f:
			subprocess.check_call(cmd.split(" "),stdout=f)

	return of

## parse options
usage = "usage: %prog [options] arg"
parser = OptionParser(usage)
parser.add_option("-a","--bedfile1",dest="bed1")
parser.add_option("-A","--bedfile2",dest="bed2")
parser.add_option("-t","--tag1",dest="tag1")
parser.add_option("-T","--tag2",dest="tag2")
parser.add_option("-b","--bamfile1",dest="bam1")
parser.add_option("-B","--bamfile2",dest="bam2")
parser.add_option("-r","--redo",dest="redo")
parser.add_option("-n","--norpmbp",dest="norpmbp")
parser.add_option("-p","--pvalue",dest="pvalue")

(options, args) = parser.parse_args()


if not options.bam2: 
	print("Usage: python3 compare_bed_files_rpmbp.py -a bedfile1 -A bedfile2 -t tag1 -T tag2 -b bamfile1 -B bamfile2 \n",file=sys.stderr)
	sys.exit()

redo=0	
if options.redo:
	redo = 1

bam=1
if options.norpmbp:
	bam=0

## read in files
## get full path for all files
files=[fpath(x) for x in [options.bed1,options.bed2,options.bam1,options.bam2]]

## move files to working directory
padd=''
if options.pvalue:
	padd='_pv'+options.pvalue

path = './rpmbp_comp_'+options.tag1+'_vs_'+options.tag2+padd
try:
    os.mkdir(path)
except OSError:
    print ("Creation of the directory %s failed" % path)
else:
    print ("Successfully created the directory %s " % path)

currentdir=os.getcwd()
os.chdir(path)
## we need to get the full pathes of each file

## make dict with sample abbreviation for files

## make sure we only take the first 5 cols of the bed files, macs1 and macs2 have different size
## if we dont do this we break the fusion script below
for i,fin in enumerate(files[0:2]):
	nof=fin+'_trunc.bed'
	print(nof)
	if redof(nof,redo):
		with open(nof,'w') as f:
			for line in open(fin):
				ls=line.strip().split()
				rest=ls[3].split("_")
				ls[3]='_'.join(rest[-2:])
				print('\t'.join(ls[0:5]),file=f)

		## set new file names to array
	files[i]=nof

myd=dict()
for i,f in enumerate(files):
	myd[i]=dict()
	myd[i]['file']=f
	if i < 2:
		myd[i]['lines']=sum(1 for line in open(files[i]))


myd[0]['ab']=options.tag1
myd[1]['ab']=options.tag2
myd[2]['ab']=options.tag1
myd[3]['ab']=options.tag2

## intersect files
k1,k2=[0,1]
of=intersectBed(k1,k2,'wo')
of1=intersectBed(k1,k2,'v')
of2=intersectBed(k2,k1,'v',1)

## fusing files now to get a real fused bed file
cmd='fuse_intervals_on_intersected_file.sh '+of
print(cmd)
run(cmd.split(" "))

ln=0
d1=dict()
d2=dict()
with open(of) as f:				
	for line in f:
		ln=ln+1
		ls=line.split()
		d1[ls[3]]=1
		d2[ls[8]]=1


f1=myd[0]['ab']
f2=myd[1]['ab']
n1=myd[0]['lines']
n2=myd[1]['lines']

statarr=[str(x) for x in [ln,0,perc(len(d1)/n1),perc(len(d2)/n2),n1,n2]]

print("{0}: {1}\t{2}: {3}\toverlap: {4}\t{5}: {6:.5f}\t{7}\t{8}: {9:.5f}\t{10}".format(f1,n1,f2,n2,ln,f1,len(d1)/n1,len(d1),f2,len(d2)/n2,len(d2) ))
print("{0} unique : {1}".format(f1,sum(1 for line in open(of1))))
print("{0} unique : {1}".format(f2,sum(1 for line in open(of2))))
print("fused regions : {0}".format(sum(1 for line in open(of+'_fused.bed'))))

## make a big bed file indicating which region is from which file, unique or intersected

if redof(of+'_tag_fused.bed',redo):
	with open (of+'_tag_fused.bed','w') as g:
		with open(of+'_fused.bed') as f:
			for line in f:
				print(line.strip()+'\tfused\t.',file=g)

		with open(of1) as f:
			for line in f:
				ls=line.strip().split()
				ls[3]=f1+'_'+ls[3]
				print('\t'.join(ls),file=g)

		with open(of2) as f:
			for line in f:
				ls=line.strip().split()
				ls[3]=f2+'_'+ls[3]
				print('\t'.join(ls),file=g)

cmd='sort -k1,1 -k2,2n '+of+'_tag_fused.bed'
bed=of+'_tag_fused_sorted.bed'

if redof(bed,redo):
	with open(bed,'w') as f:
		subprocess.call(cmd.split(),stdout=f)

rund=dict()
### now lets run the rpmbp on these files and use the cluster for it

HOME=os.getenv('HOME')
outa=[]

anyrun=0

if bam == 0:
	print("Not running rpmbp due to option B",file=stderr)
	sys.exit()

## anyways exiting here just in case
for i,bam in enumerate(files[2:4]):
	bamf=os.path.basename(bam) ## get filename only 

	## remove cluster dependence since it wont run on other clusters
	cmd='bash rpmbp.sh'
	out='rpmbp_'+bed+'_'+myd[i]['ab']
	outa.append(out)

	cmd+=' '+' '.join([bed,bam,out,'1','200'])
	print(cmd)
	cout=subprocess.run(cmd.split(),capture_output=True)

########## Plotting ##############
print("plotting now")


## now we can do the plots and then we are done
## get abbreviations
f1=options.tag1
f2=options.tag2

print(f1,f2)


wt=pd.read_csv(outa[0],sep="\t")
mt=pd.read_csv(outa[1],sep="\t")
## we do the mean on the rows already here by calling wt.mean on the Dataframe objects
df=pd.DataFrame({f1:wt.iloc[:,2],f2:mt.iloc[:,2]})
df2=df+0.0001

f1i=wt['GENE_ID'].str.contains(f1+'_') 
f2i=wt['GENE_ID'].str.contains(f2+'_')
ffi=wt['GENE_ID'].str.contains('fused_')

wt['type']=0
mt['type']=0
mt.loc[f1i,'type']=1
wt.loc[f1i,'type']=1
mt.loc[f2i,'type']=2
wt.loc[f2i,'type']=2

def plot(d,f1,f2,t1i,t2i,tfi,ti,thres,axt,hl):
	ax=sns.regplot(x=f1,y=f2,data=d.loc[t1i,:],scatter_kws={'alpha':0.15},fit_reg=False,color='red',ax=axt)
	ax=sns.regplot(x=f1,y=f2,data=d.loc[t2i,:],scatter_kws={'alpha':0.15},fit_reg=False,color='green',ax=axt)
	ax=sns.regplot(x=f1,y=f2,data=d.loc[tfi,:],scatter_kws={'alpha':0.15},fit_reg=False,color='blue',ax=axt)

	my=d.max().max().round()+1

	ax.set_xscale('log',basex=2)
	ax.set_yscale('log',basey=2)
	ax.set_xlabel(f1+'\nlog2[rpmbp+0.0001]')
	ax.set_ylabel(f2+'\nlog2[rpmbp+0.0001]')
	ax.set_ylim([thres,my])
	ax.set_xlim([thres,my])

	rho = str(df.loc[ti,:].corr(method="spearman").iloc[0,1].round(2))

	ax.set_title("union of regions\nmin rpmbp = {0}\nrho = {1}".format(thres,rho))
	ax.margins(1)
	ax.legend(handles=hl)
	return rho

r_patch = mpatches.Patch(color='red', label=f1)
g_patch = mpatches.Patch(color='green', label=f2)
b_patch = mpatches.Patch(color='blue', label='common')
hl=[r_patch, g_patch,b_patch]

sns.set()
sns.set_style("whitegrid")
fig, axs = plt.subplots(nrows=2,ncols=2,figsize=[15,15])

#filter by min rpmbp
for i,thres in enumerate([0.05,0.1,0.2,0.5]):
	ti=(df2[f1]>thres) & (df2[f2]>thres)
	t1i=(df2[f1]>thres) & (df2[f2]>thres) & f1i
	t2i=(df2[f1]>thres) & (df2[f2]>thres) & f2i
	tfi=(df2[f1]>thres) & (df2[f2]>thres) & ffi
	rho = plot(df2,f1,f2,t1i,t2i,tfi,ti,thres,axs.flatten()[i],hl)
	if thres == 0.05:
		statarr[1]=rho

plt.subplots_adjust(bottom=0.15,hspace=0.5)
#plt.savefig(f1+'_vs_'+f2+'.pdf')
plt.savefig(f1+'_vs_'+f2+'.png')

os.chdir(currentdir)
print("DONE")

harr=["Overlap","rho",",Overlap frac. w.r.t. "+f1,"Overlap frac. w.r.t. "+f2,"Peaks "+f1,"Peaks "+f2]

with open(f1+"_"+f2+"_stats.tsv",'w') as f:
	print('\t'.join(harr),file=f)
	print('\t'.join(statarr),file=f)
