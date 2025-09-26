# -----------------------------------------------------------------------------
# SPDX-FileCopyrightText: Copyright © 2025 Idiap Research Institute <contact@idiap.ch>
# SPDX-FileContributor: Raphaelle Luisier <raphaelle.luisier@idiap.ch>,Zhi Ming Xu <zhiming.xu@idiap.ch>
# SPDX-License-Identifier: GPL-3.0-only
# -----------------------------------------------------------------------------


# A. Load arguments
args       <- commandArgs(TRUE)
outdir     <- args[1]
chr        <- args[2]
region     <- args[3]
method     <- args[4]
annotation <- args[5]
limD       <- as.numeric(as.character(args[6]))
winD       <- as.numeric(as.character(args[7]))
noSeg      <- as.numeric(as.character(args[8]))
isLog      <- as.logical(as.character(args[9]))
minCov     <- as.numeric(as.character(args[10]))
tardir     <- args[11]

print(args)

print("load annotation")
load(annotation)

print("load dependencies")
source("./scripts/nested_FUN_identify_pAs.R")

library(IRanges)
library(GenomicRanges)
library(Segmentor3IsBack)

#B.Prepare annotation file
print("Prepare annotation file")
files        <- list.files(tardir)
files <- as.vector(unlist(files))
files        <- files[grep(files,pattern=".RData")]
foi          <- data.frame(files,do.call(what=rbind,args=lapply(files,function(x){
    split <- unlist(strsplit(x,split="[._]+"))
    return(split[1:(length(split)-1)])
})))
rm("files")
colnames(foi)[1]<- c("file")
colnames(foi)[ncol(foi)]<- c("chr")

# D. Import Gross Segments
print("Import gross segment")
print(region)
myregion                        <- readRDS(region)
print("region imported")
myUTR                           <- myregion[which(seqnames(myregion)==chr),]
names(myUTR)                    <- as.character(myUTR$txID)
mcols(myUTR)$transcript_id  <- myUTR$txID
out                             <- paste(outdir,"/grUTR.",chr,".RData",sep="")


#C. Load data
print("Load data")
sel         <- which(foi$chr==chr)
foi         <- foi[sel,]
myres       <-  list()
for(i in c(1:length(sel))){
    load(paste0(tardir,'/',foi$file[i]))
    myres[[i]] <- myoutput
}
metadata <- foi[,2:(ncol(foi)-1)]
names(myres)<- sapply(1:nrow(metadata),function(i) paste0(as.vector(unlist(metadata[i,,drop=T])),collapse = '.'))

#E. Identify pA sites
print("Identify pA sites")
ids      <- do.call(args=mapply(A=myres,B=names(myres),function(A,B)return(cbind(A$txID,rep(B,length(A$txID)))),SIMPLIFY=F),what=rbind)
tokeep   <- ids[,1]%in%as.character(mcols(myUTR)$transcript_id)[width(myUTR)>0]
ids      <- ids[tokeep,]
groups   <- as.factor(ids[,1])
myiso    <- list()

dir.create(outdir)
figPath <- paste0(outdir,'/Figures/')
dir.create(figPath)

logFile <- paste0(figPath,'chr',chr,'.rds')
saveRDS(list(),logFile)


for(s in c(1:length(levels(groups)))){
    print(paste("txID=",s,sep=""))
    myiso[[s]] <- tryCatch(
        {
            mergePAsTissue(tempTX=levels(groups)[s],method=method,windowsize=winD,NumSeg=noSeg,doLog=isLog,limDist=limD,MinCov=minCov,figPath = logFile)
        }
    ,error = function(cond){
            NULL
    })
}
myiso <- myiso[(!sapply(myiso,function(x) is.null(x)))]


# F. Create GRange object
myGRobject <- do.call(args=myiso,what=c)

#export.gff(object=myGRobject,con="tempGFFchr6.gtf",format="gtf")


# G. Save file
dir.create(outdir)
save(list="myGRobject",file=out)
print("Finished")
