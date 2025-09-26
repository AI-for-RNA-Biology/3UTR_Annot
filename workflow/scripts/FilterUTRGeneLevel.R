# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Raphaelle Luisier <raphaelle.luisier@idiap.ch>,Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------


library(IRanges)
library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(dplyr)

args <- commandArgs(T)
out_dir <- args[1]
myPAS_path <- args[2]
polyA_DB_path <- args[3]
polyASite_path <- args[4]
polyAModels_path <- args[5]
annot = args[6]
utr_refseq = args[7]
e75p_RefSeq = import.gff(args[8])
assembly_report = ad.table(args[9],sep = '\t'))
rpm = args[10]
APA_dir <- args[11:length(args)]

# #
# 0. Import key files
#
# 1/ Import init annotated 3'UTR
load(annot)
load(utr_refseq)
seqlevels(e75p_RefSeq) <- assembly_report$V10[match(seqlevels(e75p_RefSeq),assembly_report$V7)]

seqlevelsStyle(e75p_RefSeq) <- 'ensembl'


g3exons_ensembl <- e75p %>% filter(type=='exon' & gene_biotype %in% c('protein_coding'))
g3exons <- unique(sort(c(g3exons_ensembl,g3exons_refseq %>% filter(gene_biotype %in% c('protein_coding')))))

myPAS <- import.gff(myPAS_path,format="gtf")
polyA_DB <- import.gff(polyA_DB_path,format='gtf')
strand(polyA_DB) <- '*'
seqlevelsStyle(polyA_DB) <- 'ensembl'

polyA_Site <- import.gff(polyASite_path,format='gtf')
strand(polyA_Site) <- '*'
seqlevelsStyle(polyA_Site) <- 'ensembl'

polyA_Model <- import.gff(polyAModels_path,format='gtf')
strand(polyA_Model) <- '*'
seqlevelsStyle(polyA_Model) <- 'ensembl'


# 1. Create GRanges file for de-novo annotations
#
foi<- unlist(sapply(APA_dir,function(x) paste0(x,'/',list.files(x))))
foi <- foi[grep(foi,pattern=".RData")]
print(paste0('loading ',foi[1]))
load(foi[1])

mynewGR  <- myGRobject
start   <- start(myGRobject)
start[as.character(strand(myGRobject)) == '+'] <- start(g3utr)[match(myGRobject$txID[as.character(strand(myGRobject)) == '+'],g3utr$transcript_id)] #Correct start coordinate if it doesn't match annotation (pos strand)
end     <- end(myGRobject)
end[as.character(strand(myGRobject)) == '-'] <- end(g3utr)[match(myGRobject$txID[as.character(strand(myGRobject)) == '-'],g3utr$transcript_id)] #Correct start coordinate if it doesn't match annotation (neg strand)
chr     <- as.character(seqnames(myGRobject))
strand  <- as.character(strand(myGRobject))
isoform <- myGRobject$isoform
sampleID <- myGRobject$sampleID
transcript_id <- myGRobject$txID
gene_name <- g3utr$gene_name[match(myGRobject$txID,g3utr$transcript_id)]

for(i in c(2:length(foi))){
  print(paste0('loading ',foi[i]))
  load(foi[i])
  cur_start   <- start(myGRobject)
  cur_start[as.character(strand(myGRobject)) == '+'] <- start(g3utr)[match(myGRobject$txID[as.character(strand(myGRobject)) == '+'],g3utr$transcript_id)] #Correct start coordinate if it doesn't match annotation (pos strand)
  start  <- c(start,cur_start)
  cur_end     <- end(myGRobject)
  cur_end[as.character(strand(myGRobject)) == '-'] <- end(g3utr)[match(myGRobject$txID[as.character(strand(myGRobject)) == '-'],g3utr$transcript_id)] #Correct start coordinate if it doesn't match annotation (neg strand)
  end    <- c(end,cur_end)
  chr    <- c(chr,as.character(seqnames(myGRobject)))
  strand <- c(strand,as.character(strand(myGRobject)))
  isoform<- c(isoform,myGRobject$isoform)
  sampleID <- c(sampleID,myGRobject$sampleID)
  transcript_id <- c(transcript_id,myGRobject$txID)
  gene_name <- c(gene_name,g3utr$gene_name[match(myGRobject$txID,g3utr$transcript_id)])
  gc()
}
uniqueID <- paste0(gene_name,'.',isoform,'.',transcript_id,'.',sampleID)

mynewGR    <-  sort(GRanges(seqnames = Rle(chr),
                       ranges   = IRanges(start=start,end=end),
                       strand   = Rle(strand),
                       data.frame(uniqueID=uniqueID,sampleID=sampleID,gene_id = gene_name,gene_name=gene_name,transcript_id=transcript_id)))

#Repeat masker filter
RpM                     <- import.gff(paste0(rpm,'.gff'))
myRegion                <- mynewGR
POS                     <- as.character(strand(myRegion))=="+"
start(myRegion)[POS]    <- end(myRegion)[POS]-499
end(myRegion)[!POS]     <- start(myRegion)[!POS]+499
gOver                   <- findOverlaps(query=myRegion,subject=RpM,ignore.strand=TRUE)
idxU                    <- queryHits(gOver)
idxG                    <- subjectHits(gOver)

myStart                 <- apply(cbind(start(myRegion[idxU]),start(RpM[idxG])),1,max)
myEnd                   <- apply(cbind(end(myRegion[idxU]),end(RpM[idxG])),1,min)
mywid                   <- myEnd-myStart+1
myfracRpM               <- rep(0,length(myRegion))
myfracRpM[unique(idxU)] <- tapply(mywid,INDEX=idxU,FUN=sum)/500

torm                    <- myfracRpM>=0.3
mynewGR                    <- mynewGR[!torm,]


export.gff(object=mynewGR,con=paste0(out_dir,"hg38_Denovo.gtf"),format="gtf")
mynewGR <- import.gff(paste0(out_dir,"hg38_Denovo.gtf"))

#Add Ensembl Annotations
start   <- start(g3utr)
end     <- end(g3utr)
chr     <- as.character(seqnames(g3utr))
strand  <- as.character(strand(g3utr))
sampleID <- rep('Ensembl',length(g3utr))
transcript_id <- g3utr$transcript_id
gene_name <- g3utr$gene_name
uniqueID <- paste0(gene_name,'.',transcript_id,'.',sampleID)

ensemblGR    <-  unique(sort(GRanges(seqnames = Rle(chr),
                       ranges   = IRanges(start=start,end=end),
                       strand   = Rle(strand),
                       data.frame(uniqueID=uniqueID,sampleID=sampleID,gene_id = gene_name,gene_name=gene_name,transcript_id=transcript_id))
))


#Add RefSeq Annotations
start   <- start(g3utr_RefSeq)
end     <- end(g3utr_RefSeq)
chr     <- as.character(seqnames(g3utr_RefSeq))
strand  <- as.character(strand(g3utr_RefSeq))
sampleID <- rep('RefSeq',length(g3utr_RefSeq))
transcript_id <- g3utr_RefSeq$transcript_id
gene_name <- g3utr_RefSeq$gene
uniqueID <- paste0(gene_name,'.',transcript_id,'.',sampleID)

RefSeq_GR    <-  unique(sort(GRanges(seqnames = Rle(chr),
                       ranges   = IRanges(start=start,end=end),
                       strand   = Rle(strand),
                       data.frame(uniqueID=uniqueID,sampleID=sampleID,gene_id = gene_name,gene_name=gene_name,transcript_id=transcript_id))
))

#Merge according to order Ensembl, RefSeq, de-novo
RefSeq_Overlaps <- GenomicRanges::findOverlaps(ensemblGR,RefSeq_GR,type = 'equal')
myGR <- c(ensemblGR,RefSeq_GR[-subjectHits(RefSeq_Overlaps),])

denovo_overlaps <- GenomicRanges::findOverlaps(myGR,mynewGR,type = 'equal')
myGR <- c(myGR,mynewGR[-subjectHits(denovo_overlaps),])
myGR <- sort(myGR)

#Keep UTR with START that extend the furthest (to avoid UTR and overlaps with any other feature on the same strand)
myGR_df <- data.frame(uniqueID = myGR$uniqueID,START = start(myGR),END = end(myGR),STRAND = as.character(strand(myGR)),gene_name = myGR$gene_name,transcript_id = myGR$transcript_id) 

all_tx <- c(e75p,e75p_RefSeq)%>% dplyr::filter(type == 'transcript') %>% dplyr::select(transcript_id,tag)
all_tx_df <- data.frame(transcript_id = all_tx$transcript_id,tag = all_tx$tag)

longest_IDs <- myGR_df %>% dplyr::left_join(all_tx_df,by = c('transcript_id'='transcript_id')) %>% dplyr::filter(!(tag == 'mRNA_end_NF' | tag == 'mRNA_start_NF')) %>% group_by(gene_name) %>% filter((STRAND == '+' & START == max(START)) | (STRAND == '-' & END == min(END)) ) %>% ungroup()
myGR <- myGR %>% plyranges::filter(uniqueID %in% longest_IDs$uniqueID)

# Remove 3'UTR that extend towards any exon of another transcript on opposite strand
g3exons_flip_strand <- g3exons
strand(g3exons_flip_strand)[strand(g3exons) == "+"] <- "-"
strand(g3exons_flip_strand)[strand(g3exons) == "-"] <- "+"

utr_to_excl_exon_opp_strand <- myGR$uniqueID[unique(queryHits(findOverlaps(myGR, g3exons_flip_strand)))]
myGR <- myGR %>% plyranges::filter(!uniqueID %in% utr_to_excl_exon_opp_strand)

#add flags for apa site annotation and filter
#Focus on the regions -100;+50 around the PAS (Gruber et al, positioning)
POS                             <- as.character(strand(myGR))=="+"
subGR                           <- myGR
start(subGR)[POS]               <- end(myGR[POS])-100
end(subGR)[POS]                 <- end(myGR[POS])+50
end(subGR)[!POS]                <- start(myGR[!POS])+100
start(subGR)[!POS]              <- start(myGR[!POS])-50
start(subGR)[start(subGR)<0]    <- 1
temp.is.motif                            <- as.data.frame(queryHits(findOverlaps(query=subGR,subject=myPAS,ignore.strand=FALSE)))
is.motif                        <- rep(FALSE,length(subGR))
is.motif[unique(temp.is.motif[,1])]      <- TRUE
myGR$is.motif                     <- is.motif

temp.is.polyA_DB                            <- as.data.frame(queryHits(findOverlaps(query=subGR,subject=polyA_DB,ignore.strand=FALSE)))
is.polyA_DB                        <- rep(FALSE,length(subGR))
is.polyA_DB[unique(temp.is.polyA_DB[,1])]      <- TRUE
myGR$is.polyA_DB                     <- is.polyA_DB

temp.is.polyA_Site                            <- as.data.frame(queryHits(findOverlaps(query=subGR,subject=polyA_Site,ignore.strand=FALSE)))
is.polyA_Site                        <- rep(FALSE,length(subGR))
is.polyA_Site[unique(temp.is.polyA_Site[,1])]      <- TRUE
myGR$is.polyA_Site                     <- is.polyA_Site

temp.is.polyA_Model                            <- as.data.frame(queryHits(findOverlaps(query=subGR,subject=polyA_Model,ignore.strand=FALSE)))
is.polyA_Model                        <- rep(FALSE,length(subGR))
is.polyA_Model[unique(temp.is.polyA_Model[,1])]      <- TRUE
myGR$is.polyA_Model                     <- is.polyA_Model

export.gff(object=sort(myGR),con=paste0(out_dir,"hg38_Denovo_Ensembl_RefSeq.gtf"),format="gtf")

#Filter for motif and atleast 1 of the 3 databases. Always keep ensembl or refseq annotated
myGRFilt <- myGR %>% filter((as.logical(is.motif) == T & (as.logical(is.polyA_DB) == T | as.logical(is.polyA_Site) == T | as.logical(is.polyA_Model) == T)) | sampleID == 'Ensembl' | sampleID == 'RefSeq')

gene_incl <- data.frame(myGRFilt) %>% dplyr::group_by(gene_name) %>% dplyr::summarize(count = n()) %>% dplyr::filter(count > 1)
myGRFilt <- myGRFilt %>% filter(gene_name %in% gene_incl$gene_name)

export.gff(object=sort(myGRFilt),con=paste0(out_dir,"hg38_Denovo_Ensembl_RefSeq_Filt.gtf"),format="gtf")
export.gff(object=sort(myGRFilt) %>% select(gene_id,transcript_id = uniqueID,uniqueID=uniqueID),con=paste0(out_dir,"hg38_Denovo_Ensembl_RefSeq_Filt_Reformat.gtf"),format="gtf")

FocusEnd <- function(gr=L3,before=300,after=0){
    POS                <- as.character(strand(gr))=="+"
    grn                <- gr
    start(grn)[POS]  <- end(gr[POS])-before
    end(grn)[POS]    <- end(gr[POS])+after
    end(grn)[!POS]   <- start(gr[!POS])+before
    start(grn)[!POS] <- start(gr[!POS])-after
    start(grn)[start(grn)<0]<-1
    return(grn)
}
myGRFilt <- import.gff(paste0(out_dir,"hg38_Denovo_Ensembl_RefSeq_Filt.gtf"))
L2_300 <- FocusEnd(gr=myGRFilt,before=300,after=0)
L2_300 <- sort(L2_300) %>% plyranges::select(gene_id,transcript_id = uniqueID,uniqueID=uniqueID) 
L2_300$type <- 'exon'
export.gff(L2_300,paste0(out_dir,"hg38_Denovo_Ensembl_RefSeq_Filt_Focus_300.gtf"),format="gtf")
