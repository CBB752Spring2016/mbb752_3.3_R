GEED <- function (gctfile, clsfile, geneset) {
  
  GCT <- read.table(gctfile, sep = "\t", stringsAsFactors = FALSE, skip = 2, header = TRUE)
  GCT <- GCT[-dim(GCT)[1],] #last row is NAs --> gets rid of that row
  Names <- GCT$Name
  GCT <- GCT[,-c(1,2)] #cleaning
  
  CLS <- read.table(clsfile, skip = 2)
  
  ALL <- GCT[,CLS==0] #use CLS vector to assign to ALL
  AML <- GCT[,CLS==1] #use CLS vector to assign to AML
  pVals <- data.frame()
  for (i in (1:dim(GCT)[1])) {
    t <- t.test(ALL[i,],AML[i,])
    pVals <- rbind(pVals,t$p.value)
  }
  pVals$Name <- Names
  pVals$Name <- gsub("_[A-z]*","",pVals$Name) #formats gene names in same way as geneset.txt
  names(pVals)[1]<- "pvals"
  
  sortedP <- pVals[order(pVals$pvals),]
  
  gset <- read.table(geneset, skip = 2)
  
  sortedP$flag <- 0
  sortedP$flag[sortedP$Name %in% gset$V1] <- 1
  
  gsettotP <- sum(sortedP$pvals[sortedP$flag == 1])
  sortedP$weightedp <- 0
  sortedP$weightedp[sortedP$flag ==1] <- sortedP$pvals[sortedP$flag ==1] / gsettotP
  sortedP$phit <- 0
  sortedP$pmiss <- 0
  formulaPmiss <- 1 / (dim(sortedP)[1] - sum(sortedP$flag))
  
  for (i in (1:dim(sortedP)[1])) {
    sortedP$phit[i] <- sum(sortedP$weightedp[1:i]) #will sum all weighted p's before index
    sortedP$pmiss[i] <- formulaPmiss * i
  }
  
  sortedP$ES <- abs(sortedP$phit - sortedP$pmiss)
  
  ES <- max(sortedP$ES)
  
  return(ES)
}
  
value <- GEED("all_aml_train.gct", "all_aml_train.cls", "geneset.txt")
