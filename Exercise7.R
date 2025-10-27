library(Biostrings)
library(stringr)
library(seqinr)




DNAstring <- readDNAStringSet('seq_motif.fasta', )
motifLen <- 8
ids <- c(1,1,1,1) #??? idk fam





Score <- function(StartIdx, DNA, motifLength) {
  matrix <- c()
  
  for (i in 1:length(DNA)) {
    matrix <- c(matrix, substr(toString(DNAstring[i]), StartIdx[i],  StartIdx[i] + motifLength-1))
    
  }
  
  split_seqs <- strsplit(matrix, "")
  
  aligmatrx <- do.call(rbind, split_seqs)
  
  
  freq <- sapply(as.data.frame(aligmatrx), 
                 function(x) table(factor(x, levels = c('A', 'C', 'T', 'G'), 
                                          ordered = TRUE)))
  
  
  score_row <- apply(freq, 2, max)
  score <- 0
  for (i in 1:length(score_row)) {
    score <- score + score_row[i]
  }
  
  
  return(as.numeric(score))
  
}

#matrx <- Score(ids, DNAstring, motifLen)
#matrx








NextLeaf <- function(s, t, k) {
  
  for (i in t:1) {
    if (s[i] < k) {
      s[i] <- s[i] + 1
      return(s)
    }
    s[i] <- 1
  }
  return(s)
  
}



BFMotifSearch <- function(DNA, t, n, l) {
  s <- rep(1, length(DNA))
  if (length(s) != t) {
    return("Error 404 - too many dna sequences")
  }
  
  
  bestScore <- Score(s, DNA, l)
  while (TRUE) {
    s <- NextLeaf(s,t,n-l+1)
    if (Score(s, DNA, l) > bestScore) {
      bestScore <- Score(s,DNA,l)
      bestMotif <- s
    }
    if (identical(s, rep(1, length(DNA)))) {
      return(bestMotif)
    }
    
  }
}

bestmot <- BFMotifSearch(DNAstring, length(DNAstring), width(DNAstring)[1], 4)
print(bestmot)

