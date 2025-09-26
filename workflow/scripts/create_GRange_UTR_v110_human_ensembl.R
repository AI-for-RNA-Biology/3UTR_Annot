# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Raphaelle Luisier <raphaelle.luisier@idiap.ch>,Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------

library(GenomicFeatures)
library(Rsamtools)
library(rtracklayer)

#Set-up path
args <- commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
    config <- read_yaml('../../config/config.yml')
    e75p_path <- config$ref_gtf
    utr_path <- config$ref_annot
    utr_gtf_path <- config$utr_annot_gtf

}else{
    e75p_path <- args[[1]]
    utr_path <- args[[2]]
    utr_gtf_path <- args[[3]]
}

# A. Import Ensembl gene annotation:
e75p <- import(e75p_path)

#Keep only chromosomes and MT
e75p <- e75p[as.character(seqnames(e75p)) %in% c(seq(1,22,1),'X','Y','MT'),]


# C.1 First filter: keep only protein coding genes
F1 <- e75p$gene_biotype != "protein_coding"
# C.2 Second filter: remove when no transcript id
F2 <- is.na(e75p$transcript_id)

tokeep  <- which(!(F1|F2))
e75f    <- e75p[tokeep,]

grp             <- e75f[(strand(e75f)=="+"),]
utrp            <- grp[grp$type=="transcript",]
endUTR          <- end(ranges(grp[grp$type=="stop_codon"]))+1# if source=NMD, then there is not any STOP codon; this means I will miss all 3UTR of NMD; not sure this is a good idea.
names(endUTR)   <- as.character(mcols(grp)$transcript_id)[grp$type=="stop_codon"]
idx             <- match(as.character(mcols(utrp)$transcript_id),names(endUTR))
utrp            <- utrp[!is.na(idx)]
endUTR          <- endUTR[idx[!is.na(idx)]]
start(ranges(utrp)) <- endUTR


grn             <- e75f[(strand(e75f)=="-"),]
utrn            <- grn[grn$type=="transcript",]
endUTR          <-  start(ranges(grn[grn$type=="stop_codon"]))-1
names(endUTR)   <-  as.character(mcols(grn)$transcript_id)[grn$type=="stop_codon"]
idx             <- match(as.character(mcols(utrn)$transcript_id),names(endUTR))
utrn            <- utrn[!is.na(idx)]
endUTR          <- endUTR[idx[!is.na(idx)]]
end(ranges(utrn)) <- endUTR

utr            <- c(utrp,utrn)

 # B. Import 3'UTR using GenomicFeatures
# txdb   <- makeTxDbFromGFF(file="~/genomics/zmxu/ALS_Compartmentalization/data/GRCh38/Homo_sapiens.GRCh38.110.gtf", format="gtf")
# my3utr <- threeUTRsByTranscript(txdb, use.names=TRUE)
# mywidth<- unlist(lapply(width(my3utr),sum))

# # C. Select the end of the 3'UTR from txdb
# id1    <- as.character(mcols(utr)$transcript_id)
# id2    <- names(my3utr)
# idx    <- match(id1,id2)
# utr1   <- utr[!is.na(idx)]
# utr2   <- my3utr[idx[!is.na(idx)]]

# POS       <- as.character(strand(utr1))=="+"

# End1      <- end(ranges(utr1))
# End1[!POS]<- start(ranges(utr1))[!POS]
# End2      <- rep(NA,length(utr1))
# tempEnd   <- end(ranges(utr2))
# End2[POS] <- unlist(lapply(tempEnd[POS],max))
# tempStart <- start(ranges(utr2))
# End2[!POS]<- unlist(lapply(tempStart[!POS],min))
# diff      <- End1-End2

# Start1      <- start(ranges(utr1))
# Start1[!POS]<- end(ranges(utr1))[!POS]
# Start2      <- rep(NA,length(utr1))
# tempStart   <- start(ranges(utr2))
# Start2[POS] <- unlist(lapply(tempStart[POS],min))
# tempEnd     <- end(ranges(utr2))
# Start2[!POS]<- unlist(lapply(tempEnd[!POS],max))
# diff        <- Start1-Start2

# # D. Correct START of only those which have three nt difference
# sel                  <- width(utr1)==0
# tempStart            <- start(ranges(utr2[sel])
# Start2               <- unlist(lapply(tempStart,min))
# tempEnd              <- end(ranges(utr2[sel]))
# End2                 <- unlist(lapply(tempEnd,max))
# start(utr1)[sel]     <- apply(cbind(start(utr1[sel]),Start2),1,min)
# end(utr1)[sel]       <- apply(cbind(start(utr1[sel]),End2),1,max)
# utr[!is.na(idx)]     <- utr1

tempfun <- function(x)return(as.numeric(as.character(x)))
nUTR    <- nrow(utrn)+nrow(utrp) 
idx     <- c(match(utrp$transcript_id,e75f$transcript_id),match(utrn$transcript_id,e75f$transcript_id))

g3utr <- GRanges(seqnames =
            seqnames(e75f)[idx],
            ranges =
            IRanges(start=c(start(ranges(utrp)),start(ranges(utrn))),end=c(end(ranges(utrp)),end(ranges(utrn)))),
            strand = strand(e75f)[idx]
            )

elementMetadata(g3utr)<- elementMetadata(e75f[idx,])

gexons      <- e75f[which(e75f$type%in%c("start_codon","stop_codon","exon","UTR","CDS")),]
gtx         <- e75f[which(e75f$type%in%c("transcript")),]
gintrons    <- setdiff(gtx,gexons)
g3introns <- intersect(gintrons,g3utr)

save(list=c("e75p","e75f","g3utr",'g3introns'),file=utr_path)

g3utr$type <- 'UTR'
g3utr$gene_id <- g3utr$gene_name
export.gff(g3utr %>% filter(!is.na(gene_id)),con =utr_gtf_path,format = "gtf")
