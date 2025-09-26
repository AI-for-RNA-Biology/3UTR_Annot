# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Raphaelle Luisier <raphaelle.luisier@idiap.ch>,Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------


library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)

#Load arguments
args            <- commandArgs(TRUE)
wd              <- as.character(args[1])
chr             <- as.character(args[2])
bamfile         <- as.character(args[3])
strand          <- as.character(args[4])
pairs           <- as.character(args[5])
annotation      <- as.character(args[6])
RPM             <- as.character(args[7])
ensembl         <- as.character(args[8])

print(args)

cur_wd=getwd()

setwd(wd)
tempdir <-paste("merge_and_filter_",chr,sep="")
dir.create(tempdir)

#Merge neighbouring regions
source(paste0(cur_wd,"/scripts/mergeNeighbors.R"))

#Filter regions
source(paste0(cur_wd,"/scripts/filtering_regions.R"))

#Terminated
print("Merging and Filtering finished!")
