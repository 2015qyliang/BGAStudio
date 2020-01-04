library(tidyr)

genomeFiles <- list.files('01_GenomeFna')

filtValue <- 40

###############################################################################
if (!dir.exists('ResultSummary')) {
  dir.create('ResultSummary')
}
unlink("ResultSummary", recursive = TRUE)
dir.create('ResultSummary')

###############################################################################
###############################################################################

# COG summary
COGfuncName <- read.table('Database/COGFucName.txt', header = F, 
                          sep = '\t', quote = "",
                          stringsAsFactors = F)
colnames(COGfuncName) <- c('COG', 'class', 'function')

COGclass <- read.table('Database/COGclass.txt', header = F, 
                       sep = '\t', quote = "",
                       stringsAsFactors = F)
colnames(COGclass) <- c('class', 'ClassFunc')

mergeClassDF <- COGclass
for (fn in genomeFiles) {
  fnHeader <- strsplit(fn, '.fasta')[[1]]
  cogFilePath <- paste0('07_ResultCOG/', fnHeader, '.aa.txt')
  cogDF <- read.table(cogFilePath, header = F, 
                      sep = '\t', quote = "",
                      stringsAsFactors = F)
  cogDF <- as.vector(cogDF[which(cogDF$V3 >= filtValue),c('V2')])
  COG <- c()
  for (cogHit in cogDF) {
    cog <- rev(strsplit(cogHit, '|', fixed = T)[[1]])[1]
    COG <- append(COG, cog)
  }
  COG <- as.data.frame(COG)
  
  newCOG <- merge(COG, COGfuncName, by = 'COG', all.x = T)
  headClass <- c()
  for (cl in newCOG[,'class']) {
    CL <- strsplit(cl, '|')[[1]][1]
    headClass <- append(headClass, CL)
  }
  newCOG[,'class'] <- headClass
  
  tmpClass <- as.character(names(summary(as.factor(newCOG[,'class']))))
  tmpCount <- as.vector(summary(as.factor(newCOG[,'class'])))
  classDF <- data.frame(tmpClass, tmpCount)
  colnames(classDF) <- c('class', 'count')
  classDF <- merge(classDF, COGclass, by = 'class', all.x = T)
  colnames(classDF) <- c('class', 'count', 'function')
  
  # write.table(classDF, 
  #             paste0('ResultSummary/COG_', fnHeader, '.txt'), 
  #             quote = F,sep = '\t',
  #             row.names = F, col.names = F)
  
  mergeClassDF <- merge(mergeClassDF, classDF[, c('class', 'count')], 
                        by = 'class', all.x = T)
  colnames(mergeClassDF)[dim(mergeClassDF)[2]] <- fnHeader
  
}

numRow <- c()
for (i in c(1:dim(mergeClassDF)[1])) {
  if (sum(is.na(mergeClassDF[i,])) <= 4 ) {
    numRow <- append(numRow, i)
  }
}
mergeClassDF <- mergeClassDF[numRow, ]
for (i in 1:dim(mergeClassDF)[1]) {
  mergeClassDF[i,][is.na(mergeClassDF[i,])] <- 0
}

write.table(mergeClassDF, 
            'ResultSummary/SummaryCOG.txt', 
            quote = F,sep = '\t',
            row.names = F, col.names = T)


###############################################################################
###############################################################################

# CAZy summary
totalCAZy <- c()
for (fn in genomeFiles) {
  fnHeader <- strsplit(fn, '.fasta')[[1]]
  CAZyFilePath <- paste0('08_ResultCAZy/', fnHeader, '.aa.txt')
  CAZyDF <- read.table(CAZyFilePath, header = F, 
                       sep = '\t', quote = "",
                       stringsAsFactors = F)
  CAZyDF <- as.vector(CAZyDF[which(CAZyDF$V3 >= filtValue),c('V2')])
  CAZy <- c()
  for (cazyHit in CAZyDF) {
    cazy <- rev(strsplit(cazyHit, '|', fixed = T)[[1]])[1]
    tmpcazy <- strsplit(cazy, '|')[[1]][1]
    if (tmpcazy %in% c('A', 'C', 'G', 'P')) {
      CAZy <- append(CAZy, cazy)
    }
  }
  totalCAZy <- append(totalCAZy, CAZy)
}
CAZy <- unique(totalCAZy)

mergeCAZyDF <- data.frame(CAZy)

for (fn in genomeFiles) {
  fnHeader <- strsplit(fn, '.fasta')[[1]]
  CAZyFilePath <- paste0('08_ResultCAZy/', fnHeader, '.aa.txt')
  CAZyDF <- read.table(CAZyFilePath, header = F, 
                       sep = '\t', quote = "",
                       stringsAsFactors = F)
  CAZyDF <- as.vector(CAZyDF[which(CAZyDF$V3 >= filtValue),c('V2')])
  CAZy <- c()
  for (cazyHit in CAZyDF) {
    cazy <- rev(strsplit(cazyHit, '|', fixed = T)[[1]])[1]
    tmpcazy <- strsplit(cazy, '|')[[1]][1]
    if (tmpcazy %in% c('A', 'C', 'G', 'P')) {
      CAZy <- append(CAZy, cazy)
    }
  }
  
  tmpCAZy <- as.vector(names(summary(as.factor(CAZy))))
  tmpCount <- as.vector(summary(as.factor(CAZy)))
  CAZyDF <- data.frame(tmpCAZy, tmpCount)
  colnames(CAZyDF) <- c('CAZy', 'count')
  mergeCAZyDF <- merge(mergeCAZyDF, CAZyDF, by = 'CAZy', all.x = T)
  colnames(mergeCAZyDF)[dim(mergeCAZyDF)[2]] <- fnHeader
}

for (i in 1:dim(mergeCAZyDF)[1]) {
  mergeCAZyDF[i,][is.na(mergeCAZyDF[i,])] <- 0
}

write.table(mergeCAZyDF, 
            'ResultSummary/SummaryCAZycount.txt', 
            quote = F,sep = '\t',
            row.names = F, col.names = T)


for (i in 1:dim(mergeCAZyDF)[1]) {
  mergeCAZyDF[i,][is.na(mergeCAZyDF[i,])] <- 0
}

GHlines <- c()
GTlines <- c()
PLlines <- c()
CElines <- c()
AAlines <- c()
CBMlines <- c()
for (i in 1:dim(mergeCAZyDF)[1]) {
  cazy <- as.character(mergeCAZyDF[i,1])
  Head <- paste0(strsplit(cazy, '|')[[1]][1], strsplit(cazy, '|')[[1]][2])
  if (Head == 'GH') {
    GHlines <- append(GHlines, i)
  }
  if (Head == 'GT') {
    GTlines <- append(GTlines, i)
  }
  if (Head == 'PL') {
    PLlines <- append(PLlines, i)
  }
  if (Head == 'CE') {
    CElines <- append(CElines, i)
  }
  if (Head == 'AA') {
    AAlines <- append(AAlines, i)
  }
  if (Head == 'CB') {
    CBMlines <- append(CBMlines, i)
  }
}

GHs <- colSums(mergeCAZyDF[GHlines,c(2:6)])
GTs <- colSums(mergeCAZyDF[GTlines,c(2:6)])
PLs <- colSums(mergeCAZyDF[PLlines,c(2:6)])
CEs <- colSums(mergeCAZyDF[CElines,c(2:6)])
AAs <- colSums(mergeCAZyDF[AAlines,c(2:6)])
CBMs <- colSums(mergeCAZyDF[CBMlines,c(2:6)])

sumColsDF <- data.frame(GHs, GTs, PLs, CEs, AAs, CBMs)
rownames(sumColsDF)
colnames(sumColsDF)

write.table(sumColsDF, 
            'ResultSummary/SummaryCAZycountSums.txt', 
            quote = F,sep = '\t',
            row.names = T, col.names = T)


###############################################################################
###############################################################################

# OG5 summary
totalOG5 <- c()
for (fn in genomeFiles) {
  fnHeader <- strsplit(fn, '.fasta')[[1]]
  og5FilePath <- paste0('06_ResultOrthoMclOG5/', fnHeader, '.aa.txt')
  og5DF <- read.table(og5FilePath, header = F, 
                      sep = '\t', quote = "",
                      stringsAsFactors = F)
  og5DF <- as.vector(og5DF[which(og5DF$V3 >= filtValue),c('V2')])
  og5 <- c()
  for (ogHit in og5DF) {
    og <- rev(strsplit(ogHit, '|', fixed = T)[[1]])[1]
    if (og != 'no_group') {
      og5 <- append(og5, og)
    }
  }
  totalOG5 <- append(totalOG5, og5)
}
OG5 <- unique(totalOG5)

mergeOG5DF <- data.frame(OG5)

for (fn in genomeFiles) {
  fnHeader <- strsplit(fn, '.fasta')[[1]]
  og5FilePath <- paste0('06_ResultOrthoMclOG5/', fnHeader, '.aa.txt')
  og5DF <- read.table(og5FilePath, header = F, 
                      sep = '\t', quote = "",
                      stringsAsFactors = F)
  og5DF <- as.vector(og5DF[which(og5DF$V3 >= filtValue),c('V2')])
  og5 <- c()
  for (ogHit in og5DF) {
    og <- rev(strsplit(ogHit, '|', fixed = T)[[1]])[1]
    if (og != 'no_group') {
      og5 <- append(og5, og)
    }
  }
  
  tmpOG5 <- as.vector(names(summary(as.factor(og5))))
  tmpCount <- as.vector(summary(as.factor(og5)))
  ogDF <- data.frame(tmpOG5, tmpCount)
  colnames(ogDF) <- c('OG5', 'count')
  mergeOG5DF <- merge(mergeOG5DF, ogDF, by = 'OG5', all.x = T)
  colnames(mergeOG5DF)[dim(mergeOG5DF)[2]] <- fnHeader
}

numRow <- c()
for (i in c(1:dim(mergeOG5DF)[1])) {
  if (sum(is.na(mergeOG5DF[i,])) <= 4 ) {
    numRow <- append(numRow, i)
  }
}
mergeOG5DF <- mergeOG5DF[numRow, ]

for (i in 1:dim(mergeOG5DF)[1]) {
  mergeOG5DF[i,][is.na(mergeOG5DF[i,])] <- 0
}

write.table(mergeOG5DF, 
            'ResultSummary/SummaryOG5count.txt', 
            quote = F,sep = '\t',
            row.names = F, col.names = T)


###############################################################################
###############################################################################



