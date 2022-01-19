# Trim sequences and delete short sequences.
L <- as.numeric(readline("\nPlease enter the trimming sequence length (e.g. 250):\t"))

feature.DaDa2 <-row.names(seq.DaDa2)
seq.DaDa2 <- cbind(feature.DaDa2, seq.DaDa2)
colnames(seq.DaDa2)[2] <- "sequences"
length <- nchar(seq.DaDa2$sequences)
seq.DaDa2 <- cbind(seq.DaDa2, length)
colnames(seq.DaDa2)[3] <- "length"
seq.DaDa2 <- subset(seq.DaDa2, seq.DaDa2$length >= L)
seq.DaDa2$sequences <- substring(seq.DaDa2$sequences, 1, L)
seq.DaDa2$length.trim <- nchar(seq.DaDa2$sequences)
seq.DaDa2 <- seq.DaDa2[, -3] 
colnames(seq.DaDa2)[3] <- "length" 

feature.Deblur <-row.names(seq.Deblur)
seq.Deblur <- cbind(feature.Deblur, seq.Deblur)
colnames(seq.Deblur)[2] <- "sequences"
length <- nchar(seq.Deblur$sequences)
seq.Deblur <- cbind(seq.Deblur, length)
colnames(seq.Deblur)[3] <- "length"
seq.Deblur$sequences <- substring(seq.Deblur$sequences, 1, L)
seq.Deblur$length.trim <- nchar(seq.Deblur$sequences)
seq.Deblur <- seq.Deblur[, -3]
colnames(seq.Deblur)[3] <- "length" 

feature.UNOISE3 <-row.names(seq.UNOISE3)
seq.UNOISE3 <- cbind(feature.UNOISE3, seq.UNOISE3)
colnames(seq.UNOISE3)[2] <- "sequences"
length <- nchar(seq.UNOISE3$sequences)
seq.UNOISE3 <- cbind(seq.UNOISE3, length)
colnames(seq.UNOISE3)[3] <- "length"
seq.UNOISE3$sequences <- substring(seq.UNOISE3$sequences, 1, L)
seq.UNOISE3$length.trim <- nchar(seq.UNOISE3$sequences)
seq.UNOISE3 <- seq.UNOISE3[, -3]
colnames(seq.UNOISE3)[3] <- "length" 

feature.UPARSE <-row.names(seq.UPARSE)
seq.UPARSE <- cbind(feature.UPARSE, seq.UPARSE)
colnames(seq.UPARSE)[2] <- "sequences"
length <- nchar(seq.UPARSE$sequences)
seq.UPARSE <- cbind(seq.UPARSE, length)
colnames(seq.UPARSE)[3] <- "length"
seq.UPARSE$sequences <- substring(seq.UPARSE$sequences, 1, L)
seq.UPARSE$length.trim <- nchar(seq.UPARSE$sequences)
seq.UPARSE <- seq.UPARSE[, -3]
colnames(seq.UPARSE )[3] <- "length" 

# Output data.
dir.create("sequences")

write.table(seq.DaDa2,file ="sequences/seq.DaDa2.txt")
write.table(seq.Deblur,file ="sequences/seq.Deblur.txt")
write.table(seq.UNOISE3,file ="sequences/seq.UNOISE3.txt")
write.table(seq.UPARSE,file ="sequences/seq.UPARSE.txt")
