# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Raphaelle Luisier <raphaelle.luisier@idiap.ch>,Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------


#Merge all GTF files in a directory
library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

args <- commandArgs(T)
input_dir <- args[1]
out_file <- args[2]

foi <- list.files(input_dir)
foi <- foi[grep(foi,pattern="A.gtf")]
myPAS <- do.call(lapply(foi,function(x)return(import.gff(con=paste(input_dir,x,sep=""),format="gtf"))),what=c)
export.gff(object=myPAS,con=out_file,format="gtf")
