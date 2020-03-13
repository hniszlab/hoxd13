#!/usr/bin/Rscript

# generate_figure_7B.R -  R script
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

DPATH="../../Fdata/"

load(paste(DPATH,"Figure7/pmat.tf2_10.RData",sep=""))
load(paste(DPATH,"Figure7/pmat.dbd_10.RData",sep=""))


### Figure7B - pondr
svg("Figure7B.svg",3,6)
m=list(pmat.dbd_10$matrix[,1],pmat.tf2_10$matrix[,1])
names(m)=c("DBDs","IDRs")
boxplot(m,col=c("grey","blue"),las=1,ylab="PONDR score")
dev.off()

write.table(pmat.dbd_10$matrix[,1],"PONDR_scores_DBD.tsv",col.names=F,quote=F)
write.table(pmat.tf2_10$matrix[,1],"PONDR_scores_TF.tsv",col.names=F,quote=F)

