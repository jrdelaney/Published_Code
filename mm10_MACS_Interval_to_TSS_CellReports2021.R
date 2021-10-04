library(dplyr)
library(data.table)

#set distance from TSS (in bp)
priorTSS <- 250
postTSS <- 1000

#mm10 data from UCSC table browser RefGene
mm10genes <- as.data.frame(fread(paste0("Annotations/mm10refGene.txt"), header=TRUE, skip=0, stringsAsFactors=FALSE, sep="\t"))

mm10positive <- subset(mm10genes, strand=="+")[,c("chrom", "txStart","txEnd","strand","name2")]
mm10negative <- subset(mm10genes, strand=="-")[,c("chrom", "txStart","txEnd","strand","name2")]

mm10positive$lowRange <- mm10positive$txStart - priorTSS
mm10positive$highRange <- mm10positive$txStart +postTSS
mm10negative$lowRange <- mm10negative$txEnd - postTSS
mm10negative$highRange <- mm10negative$txEnd + priorTSS

mm10positive <- mm10positive %>%
  rowwise() %>%
  mutate(chromExtra = nchar(chrom))
mm10positive <- mm10positive %>%
  filter(chromExtra<6)
mm10positive <- mm10positive[,c("chrom", "txStart","txEnd","strand","name2","lowRange","highRange")]

mm10negative <- mm10negative %>%
  rowwise() %>%
  mutate(chromExtra = nchar(chrom))
mm10negative <- mm10negative %>%
  filter(chromExtra<6)
mm10negative <- mm10negative[,c("chrom", "txStart","txEnd","strand","name2","lowRange","highRange")]

mm10ranges<- rbind(mm10positive, mm10negative)
#rm(mm10negative,mm10positive,mm10genes)

backgroundGeneList <- unique(mm10ranges$name2)
write.table(backgroundGeneList, "ChIPseq/ChIP Output/geneBackgroundList.csv", row.names = FALSE, col.names = FALSE,quote = FALSE, sep = ",")

#read interval files from MACS outputs
intervalFiles = list.files("ChIPseq/Intervals/", pattern="*.interval")
intervalFileLocations <- as.list(rep("ChIPseq/Intervals/",length(intervalFiles)))
for (file in 1:length(intervalFiles)){
  intervalFileLocations[file] <- paste(intervalFileLocations[file],intervalFiles[file], sep="")
}

intervalLabels <- gsub(".interval", "",intervalFiles)

totalPeaks <- 0
intervalTables <- list()
for(i in 1:length(intervalFiles)){
  intervalTables[[i]] <- as.data.frame(fread(intervalFileLocations[[i]], header=TRUE, skip=18))
totalPeaks <- totalPeaks + nrow(intervalTables[[i]])
  }

mm10ranges$aboveLow <- rep(FALSE,nrow(mm10ranges))
mm10ranges$belowHigh <- rep(FALSE,nrow(mm10ranges))
mm10ranges$both <- rep(FALSE,nrow(mm10ranges))

chromNum <- "chr"
gene <- NULL
geneHitList <- NULL
betweenCk <- mm10ranges
peakNum <- 0
geneList <- rep(list(NULL),length(intervalFiles))
score_table <- data.frame(name2 = "", score = 0)
score_table_list <- NULL
score_table_byInterval <- NULL
score_table_byInterval_list <- NULL

for(intervalFile in 1:length(intervalFiles)){ #intervalFile <- 1
geneHitList <- NULL
one <- intervalTables[[intervalFile]]
one$peak <- intervalTables[[intervalFile]]$start + intervalTables[[intervalFile]]$summit
colnames(one) <- c('chrom',colnames(one)[2:length(colnames(one))])
score_table <- data.frame(name2 = "", score = 0)
score_table_list <- NULL
score_table_byInterval <- NULL


numRowChrom <- 0
for(chromosome in 1:length(unique(one$chrom))){ #chromosome <- 1
  oneData <- subset(one,chrom== unique(one$chrom)[chromosome])
  chromNum <- oneData[1,1]
  numRowChrom <- nrow(oneData)

for(peak in 1:numRowChrom){ #peak <- 1
  if(peak==1){chromNum <- oneData[peak,1]}
  betweenCk <- mm10ranges[which(mm10ranges$chrom==chromNum),]
  betweenCk$both <- (betweenCk$lowRange < oneData$peak[peak]) & (betweenCk$highRange > oneData$peak[peak]) | #checks if peak is within set range
    (betweenCk$lowRange > oneData$start[peak]) & (betweenCk$highRange < oneData$end[peak]) | #or if entire peak spans range
    (betweenCk$lowRange < oneData$start[peak]) & (betweenCk$highRange > oneData$start[peak]) | #or if peak edge (start) in range
    (betweenCk$lowRange < oneData$end[peak]) & (betweenCk$highRange > oneData$end[peak]) #or if peak edge (end) in range
  betweenCk$everBoth <- betweenCk$lowRange < -10000000000
  betweenCk$everBoth <- betweenCk$everBoth | betweenCk$both
  if(peak==numRowChrom){geneHitList <- c(geneHitList,unique(unlist(betweenCk[betweenCk$everBoth,"name2"])))}
  if(length(unique(unlist(betweenCk[betweenCk$everBoth,"name2"]))) > 0){
  score_table <- data.frame(name2 = unique(unlist(betweenCk[betweenCk$everBoth,"name2"]))
                            , score = oneData$`-10*log10(pvalue)`[peak])
  score_table_list[[length(score_table_list)+1]] <- score_table
  }
  }
}

geneList[[intervalFile]] <- unique(geneHitList)
score_table_byInterval <- do.call("rbind",score_table_list)
write.table( score_table_byInterval
            ,paste0("ChIPseq/ChIP Output/",intervalLabels[intervalFile]," geneScores.csv")
            ,row.names = FALSE,col.names=TRUE, sep=",", quote = FALSE)
if(intervalFile == 1){
score_table_byInterval_list[[1]] <- score_table_byInterval
} else {
  score_table_byInterval_list[[intervalFile]] <- score_table_byInterval
}

}

geneHitList <- unique(unlist(
  lapply(1:length(intervalFiles), function(intervalFile){ #intervalFile <- 1
    return(unique(score_table_byInterval_list[[intervalFile]][,"name2"]))
  })
))
write.table(geneHitList, "ChIPseq/ChIP Output/geneHitList.csv", row.names = FALSE,col.names=FALSE, sep=",", quote = FALSE)

geneScoreMeans <- 
  matrix(0, nrow = length(geneHitList), ncol = length(intervalFiles))

row.names(geneScoreMeans) <- geneHitList
colnames(geneScoreMeans) <- intervalLabels

for(gene in geneHitList){ # gene <- geneHitList[1]
  geneScoreMeans[which(geneHitList == gene),] <- unlist(
  lapply(1:length(intervalFiles), function(intervalFile){ #intervalFile <- 1
    score_table_byInterval <- score_table_byInterval_list[[intervalFile]]
    return(
      max(
      score_table_byInterval[which(score_table_byInterval$name2 == gene),"score"]
      ,0,na.rm = TRUE))
  })
  )
}

geneScoreMeans[which(geneScoreMeans == 0)] <- NA
geneScoreMeans_mean <- rowMeans(geneScoreMeans, na.rm = TRUE)
geneScoreMeans <- as.data.frame(geneScoreMeans)
geneScoreMeans$mean <- geneScoreMeans_mean
geneScoreMeans$gene_symbol <- row.names(geneScoreMeans)
geneScoreMeans <- geneScoreMeans[,c("gene_symbol","mean",intervalLabels)]
write.table(geneScoreMeans, "ChIPseq/ChIP Output/geneScores.csv", row.names = FALSE,col.names=TRUE, sep=",", quote = FALSE)

