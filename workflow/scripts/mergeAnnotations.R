# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Raphaelle Luisier <raphaelle.luisier@idiap.ch>,Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------


library(grDevices)
library(Rsamtools)
library(IRanges)
library(GenomicRanges)
library(rtracklayer)

args <- commandArgs(T)


out <- args[1]
annot <- args[2]

# # A. Import new annotation
foi1_foi2 <- args[3:length(args)]
foi1    <- foi1_foi2[sapply(foi1_foi2,function(x) grepl(pattern = '.6.gff',x=x))]
foi2    <- foi1_foi2[sapply(foi1_foi2,function(x) grepl(pattern = '.extended.gff',x=x))]
fois    <- c(foi1,foi2)

myS     <- do.call(lapply(fois,function(x)return(import.gff(con=x,version = '1'))),what=c)

# B. Import known 3' UTR and annotation of transcript
load(annot)
#group variable not in annot, need to export then import (don't know why...)
export.gff(g3utr,con = paste0(out,'.g3utr.gff'),format = 'gff3')
g3utr_v1 <- import.gff(paste0(out,'.g3utr.gff'),version = '1')
system(paste0('rm ',out,'.g3utr.gff'))

# C. Extract the longest 3'UTR isoforms from new annotation
ids          <- as.character(mcols(g3utr)$transcript_id)
group        <- as.character(mcols(g3utr_v1)$group)
txID         <- ids[match(as.character(mcols(myS)$group),group)]
myS$txID     <- txID
dL           <- width(myS)-width(g3utr)[match(myS$txID,ids)]

newS         <- GRanges(seqnames = c(seqnames(myS),seqnames(g3utr)),
                       ranges   = c(ranges(myS),ranges(g3utr)),
                       strand   = c(strand(myS),strand(g3utr)),
                       txID     = c(myS$txID,as.character(mcols(g3utr)$transcript_id)),
                       dL       = c(dL,rep(0,length(g3utr)))
           )

no.iso       <- as.vector(table(newS$txID))[match(newS$txID,names(table(newS$txID)))]
newS$no.iso  <- no.iso


sub1        <- newS[newS$no.iso==1,]
sub2        <- newS[newS$no.iso>1,]
myTX        <- unique(sub2$txID)
temp        <- cbind(sub2$dL,sub2$txID,c(1:length(sub2)))
theLongest  <- do.call(lapply(myTX, function(x){
                            ix  <- which(temp[,2]==x)
                            sub <- temp[ix,]
                            return(sub[sort(as.numeric(sub[,1]),decreasing=T,index.return=T)$ix[1],3])
                         }
                         ),what=c)

myLongest              <- c(sub1,sub2[as.numeric(theLongest),])
POS                    <- as.character(strand(myLongest))=="+"
end(myLongest)[POS]    <- end(myLongest)[POS] + 100
start(myLongest)[!POS] <- start(myLongest)[!POS] - 100
start(myLongest)[start(myLongest)<1]<-1

export.gff(object=myLongest,con=out,format="gtf")

#Create bed file
start        <- start(myLongest)
end          <- end(myLongest)
chr          <- as.character(seqnames(myLongest))
name         <- myLongest$txID
strand       <- as.character(strand(myLongest))
outUTR          <- cbind(chr, start,end,name,score=rep(1,length(myLongest)),strand)
write.table(x=outUTR,file=gsub(x=out,pattern = 'gtf',replacement = 'bed'), append = FALSE, quote = FALSE, sep = "\t", eol = "\n", na = "NA", row.names = FALSE,col.names = FALSE)

#Create RDS file
saveRDS(myLongest,gsub(x=out,pattern = '.gtf',replacement = '.rds'))