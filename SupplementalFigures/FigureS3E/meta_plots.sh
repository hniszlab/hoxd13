#!/bin/bash

# make_plots.sh -  Shell script
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

U='Usage: meta_plots.sh   bedfile  bamfile NAME [options] flanksize   bins   flankregionsize'

if [[ -z $3 ]];then
	echo 
	echo No NAME given as argv[3] 
	echo 
	echo $U
	echo bedfile          bedfile to get rpmbp signal for
	echo bamfile          bamfile to calculate signal for
	echo NAME             name of output folder
	echo 
	echo [options]
	echo flanksize        size of bedregion extension e.g. 1000
	echo bin              default: 50, number of bins a signal is gotten for. It is mainly for plotting purposes
	echo flankregionsize  default: 1000, size of the region left and right of the bed region to compare with
	echo 
	exit
fi

## path to BradnerPipeline
P=~/a/tools/pipeline

BEDFILE=$(readlink -f $1)
BAM=$(readlink -f $2)
NAME=$3
FLANK=1000
BINS=50
FLANKREGION=3000
let BS=$FLANKREGION/$BINS

if [[ $4 ]];then
	FLANK=$4
fi

if [[ $5 ]];then 
	BINS=$5
fi

if [[ $6 ]];then
	FLANKREGION=$6
fi

LOG=run.log

## create output directory
bn=$(basename $BEDFILE .bed)
OUTD=rpmbp_$NAME
if [[ ! -d $OUTD ]];then
	mkdir $OUTD
fi
cd $OUTD

BIN5ADD=10 ## bins to add left and right of a region. Cut off later when plotting again.
let BIN5ADD2=2*$BIN5ADD

## make a gff file with expanded regions
python3 ../bed2gff.py $BEDFILE $FLANK > ${bn}_${FLANK}.gff

## here we add the extra bins which we cut off later again
BINS=$BINS BINA=$BIN5ADD perl -ane '$bins=$ENV{BINS};$bs=int(($F[4]-$F[3]+1)/$bins);print "$F[0]\t$F[1]\t$F[2]\t",$F[3]-($bs*$ENV{BINA}),"\t",$F[4]+($bs*$ENV{BINA}),"\t$F[5]\t$F[6]\t$F[7]\t$F[8];binsize=$bs;orig=${F[0]}_${F[3]}_${F[4]}\n";' ${bn}_${FLANK}.gff >  ${bn}_center_expanded.gff

## gff3 has all coords included so we add or subtract 1 depending on the side
## add 300nt to the flanks for smooting of the plots later, this is 6 bins of size 50 bp
BINSFS=$BINS

FL=$FLANKREGION NB=$BINS BINA=$BIN5ADD perl -ane '$fl=$ENV{"FL"};$badd=$ENV{BINA}*$fl/$ENV{NB};$F[4]=($F[3]-1)+$badd;$F[3]=$F[3]-$fl-$badd;$F[2]="flankL";print join("\t",@F),"\n";' ${bn}_${FLANK}.gff > ${bn}_left_${FLANKREGION}_flank.gff
FL=$FLANKREGION NB=$BINS BINA=$BIN5ADD perl -ane '$fl=$ENV{"FL"};$badd=$ENV{BINA}*$fl/$ENV{NB};$F[3]=($F[4]+1)-$badd;$F[4]=$F[4]+$fl+$badd;$F[2]="flankR";print join("\t",@F),"\n";' ${bn}_${FLANK}.gff > ${bn}_right_${FLANKREGION}_flank.gff

let BINSF=$BINS+$BIN5ADD2
let BINSC=$BINS+$BIN5ADD2 ## five at each side which we remove later on again 

## now run rpmbp.sh on it
## once it is finished run the plot_making
## option -e extends reads by 200 bp on each side before mapping it to a bin
EXT=200 ## this option can usually be kept
#

if [[ ! -f "success" ]];then 
	echo $P/bamToGFF_turbo.py -b $BAM -i ${bn}_center_expanded.gff -o rpmbp_${BINS}bins_${EXT}ext_center_${bn} -s both -e $EXT -m ${BINSC} -r \& >> run.log
	$P/bamToGFF_turbo.py -b $BAM -i ${bn}_center_expanded.gff -o rpmbp_${BINS}bins_${EXT}ext_center_${bn} -s both -e $EXT -m ${BINSC} -r &
	PID1=$!

	echo $P/bamToGFF_turbo.py -b $BAM -i ${bn}_left_${FLANKREGION}_flank.gff  -o rpmbp_${BINS}bins_${EXT}ext_left_${bn}   -s both -e $EXT -m ${BINSF} -r \& >> run.log
	$P/bamToGFF_turbo.py -b $BAM -i ${bn}_left_${FLANKREGION}_flank.gff  -o rpmbp_${BINS}bins_${EXT}ext_left_${bn}   -s both -e $EXT -m ${BINSF} -r &
	PID2=$!

	echo $P/bamToGFF_turbo.py -b $BAM -i ${bn}_right_${FLANKREGION}_flank.gff -o rpmbp_${BINS}bins_${EXT}ext_right_${bn}  -s both -e $EXT -m ${BINSF} -r \& >> run.log
	$P/bamToGFF_turbo.py -b $BAM -i ${bn}_right_${FLANKREGION}_flank.gff -o rpmbp_${BINS}bins_${EXT}ext_right_${bn}  -s both -e $EXT -m ${BINSF} -r &
	PID3=$!

	## so wait until all jobs finished
	echo perl ../check_running_process.pl $PID1 $PID2 $PID3 >> run.log
	perl ../check_running_process.pl $PID1 $PID2 $PID3
fi

touch success

## now parse the files to tsv and read into R to make the plots
python3 ../sanity_check.py rpmbp_${BINS}bins_${EXT}ext_center_${bn}
python3 ../sanity_check.py rpmbp_${BINS}bins_${EXT}ext_left_${bn}
python3 ../sanity_check.py rpmbp_${BINS}bins_${EXT}ext_right_${bn}

## now plot rpmbp in R
echo Rscript make_plot.R --rpmbp=rpmbp_${BINS}bins_${EXT}ext_center_${bn}.tsv --name=$bn --bins=$BIN5ADD --flank=$FLANKREGION >> run.log
Rscript ../make_plot.R --rpmbp=rpmbp_${BINS}bins_${EXT}ext_center_${bn}.tsv --name=$bn --bins=$BIN5ADD --flank=$FLANKREGION --title=${NAME}

echo 
echo DONE
echo
cd ..
