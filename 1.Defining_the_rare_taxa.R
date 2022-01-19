# Calculate relative abundance.
abunda.DaDa2 <- table.DaDa2 / colSums(table.DaDa2)
abunda.Deblur <- table.Deblur / colSums(table.Deblur)
abunda.UNOISE3 <- table.UNOISE3 / colSums(table.UNOISE3)
abunda.UPARSE <- table.UPARSE / colSums(table.UPARSE)

# Delete the feature that occurs once in triplicate.
Occurrence <- function(x) sum(x != 0)

m <- as.numeric(readline("\nPlease enter the number of the samples (e.g. 3):\t"))
n <- as.numeric(readline("\nPlease enter the occurrence of the features among all samples (e.g. 2):\t"))

abunda.DaDa2$Occurrence <- apply(abunda.DaDa2, 1, Occurrence)
pick.DaDa2 <- subset(abunda.DaDa2, abunda.DaDa2$Occurrence >= n)
exclude.DaDa2 <- subset(abunda.DaDa2, abunda.DaDa2$Occurrence < n)

abunda.Deblur$Occurrence <- apply(abunda.Deblur, 1, Occurrence)
pick.Deblur <- subset(abunda.Deblur, abunda.Deblur$Occurrence >= n)
exclude.Deblur <- subset(abunda.Deblur, abunda.Deblur$Occurrence < n)

abunda.UNOISE3$Occurrence <- apply(abunda.UNOISE3, 1, Occurrence)
pick.UNOISE3 <- subset(abunda.UNOISE3, abunda.UNOISE3$Occurrence >= n)
exclude.UNOISE3 <- subset(abunda.UNOISE3, abunda.UNOISE3$Occurrence < n)

abunda.UPARSE$Occurrence <- apply(abunda.UPARSE, 1, Occurrence)
pick.UPARSE <- subset(abunda.UPARSE, abunda.UPARSE$Occurrence >= n)
exclude.UPARSE <- subset(abunda.UPARSE, abunda.UPARSE$Occurrence < n)

pick.DaDa2 <- subset(pick.DaDa2,select=-Occurrence)
pick.Deblur <- subset(pick.Deblur,select=-Occurrence)
pick.UNOISE3 <- subset(pick.UNOISE3,select=-Occurrence)
pick.UPARSE <- subset(pick.UPARSE,select=-Occurrence)

# Average the abundance of each feature.
pick.DaDa2$mean.abundance  <- apply(pick.DaDa2, 1, mean)
pick.Deblur$mean.abundance  <- apply(pick.Deblur, 1, mean)
pick.UNOISE3$mean.abundance  <- apply(pick.UNOISE3, 1, mean)
pick.UPARSE$mean.abundance  <- apply(pick.UPARSE, 1, mean)

# Define the abundant (>1%), transitional (0.1~1%) and rare (<0.1%) taxa.
pick.DaDa2 <- as.data.frame(pick.DaDa2[order(pick.DaDa2$mean.abundance),])
pick.Deblur <- as.data.frame(pick.Deblur[order(pick.Deblur$mean.abundance),])
pick.UNOISE3 <- as.data.frame(pick.UNOISE3[order(pick.UNOISE3$mean.abundance),])
pick.UPARSE <- as.data.frame(pick.UPARSE[order(pick.UPARSE$mean.abundance),])

abundant.DaDa2 <- pick.DaDa2[which(pick.DaDa2$mean.abundance > 0.01),]
transitional.DaDa2 <- pick.DaDa2[which(pick.DaDa2$mean.abundance >= 0.001 & pick.DaDa2$mean.abundance <= 0.01),]
rare.DaDa2 <- pick.DaDa2[which(pick.DaDa2$mean.abundance < 0.001),]

abundant.Deblur <- pick.Deblur[which(pick.Deblur$mean.abundance > 0.01),]
transitional.Deblur <- pick.Deblur[which(pick.Deblur$mean.abundance >= 0.001 & pick.Deblur$mean.abundance <= 0.01),]
rare.Deblur <- pick.Deblur[which(pick.Deblur$mean.abundance < 0.001),]

abundant.UNOISE3 <- pick.UNOISE3[which(pick.UNOISE3$mean.abundance > 0.01),]
transitional.UNOISE3 <- pick.UNOISE3[which(pick.UNOISE3$mean.abundance >= 0.001 & pick.UNOISE3$mean.abundance <= 0.01),]
rare.UNOISE3 <- pick.UNOISE3[which(pick.UNOISE3$mean.abundance < 0.001),]

abundant.UPARSE <- pick.UPARSE[which(pick.UPARSE$mean.abundance > 0.01),]
transitional.UPARSE <- pick.UPARSE[which(pick.UPARSE$mean.abundance >= 0.001 & pick.UPARSE$mean.abundance <= 0.01),]
rare.UPARSE <- pick.UPARSE[which(pick.UPARSE$mean.abundance < 0.001),]

# Calculate the number and abundance of each taxa.
abundant.DaDa2$feature <- rownames(abundant.DaDa2)
transitional.DaDa2$feature <- rownames(transitional.DaDa2)
rare.DaDa2$feature <- rownames(rare.DaDa2)

abundant.Deblur$feature <- rownames(abundant.Deblur)
transitional.Deblur$feature <- rownames(transitional.Deblur)
rare.Deblur$feature <- rownames(rare.Deblur)

abundant.UNOISE3$feature <- rownames(abundant.UNOISE3)
transitional.UNOISE3$feature <- rownames(transitional.UNOISE3)
rare.UNOISE3$feature <- rownames(rare.UNOISE3)

abundant.UPARSE$feature <- rownames(abundant.UPARSE)
transitional.UPARSE$feature <- rownames(transitional.UPARSE)
rare.UPARSE$feature <- rownames(rare.UPARSE)

table.DaDa2$feature <- rownames(table.DaDa2)
table.Deblur$feature <- rownames(table.Deblur)
table.UNOISE3$feature <- rownames(table.UNOISE3)
table.UPARSE$feature <- rownames(table.UPARSE)

reads.abundant.DaDa2  <-merge(abundant.DaDa2, table.DaDa2, by.x = "feature", by.y = "feature") 
reads.transitional.DaDa2 <-merge(transitional.DaDa2 , table.DaDa2, by.x = "feature", by.y = "feature") 
reads.rare.DaDa2 <-merge(rare.DaDa2 ,table.DaDa2, by.x = "feature", by.y = "feature") 

reads.abundant.Deblur  <-merge(abundant.Deblur, table.Deblur, by.x = "feature", by.y = "feature") 
reads.transitional.Deblur <-merge(transitional.Deblur , table.Deblur, by.x = "feature", by.y = "feature") 
reads.rare.Deblur <-merge(rare.Deblur ,table.Deblur, by.x = "feature", by.y = "feature") 

reads.abundant.UNOISE3  <-merge(abundant.UNOISE3, table.UNOISE3, by.x = "feature", by.y = "feature") 
reads.transitional.UNOISE3 <-merge(transitional.UNOISE3 , table.UNOISE3, by.x = "feature", by.y = "feature") 
reads.rare.UNOISE3 <-merge(rare.UNOISE3 ,table.UNOISE3, by.x = "feature", by.y = "feature") 

reads.abundant.UPARSE  <-merge(abundant.UPARSE, table.UPARSE, by.x = "feature", by.y = "feature") 
reads.transitional.UPARSE <-merge(transitional.UPARSE , table.UPARSE, by.x = "feature", by.y = "feature") 
reads.rare.UPARSE <-merge(rare.UPARSE ,table.UPARSE, by.x = "feature", by.y = "feature") 

sum.abundant.DaDa2 <- colSums(subset(reads.abundant.DaDa2,select=-c(feature, mean.abundance)))
sum.transitional.DaDa2 <- colSums(subset(reads.transitional.DaDa2,select=-c(feature, mean.abundance)))
sum.rare.DaDa2 <- colSums(subset(reads.rare.DaDa2,select=-c(feature, mean.abundance)))
sum.DaDa2 <- t(cbind.data.frame(sum.abundant.DaDa2, sum.transitional.DaDa2, sum.rare.DaDa2))
sum.DaDa2 <- rbind.data.frame(sum.DaDa2, colSums(sum.DaDa2))
row.names(sum.DaDa2)[4] <- "sum"

sum.abundant.Deblur <- colSums(subset(reads.abundant.Deblur,select=-c(feature, mean.abundance)))
sum.transitional.Deblur <- colSums(subset(reads.transitional.Deblur,select=-c(feature, mean.abundance)))
sum.rare.Deblur <- colSums(subset(reads.rare.Deblur,select=-c(feature, mean.abundance)))
sum.Deblur <- t(cbind.data.frame(sum.abundant.Deblur, sum.transitional.Deblur, sum.rare.Deblur))
sum.Deblur <- rbind.data.frame(sum.Deblur, colSums(sum.Deblur))
row.names(sum.Deblur)[4] <- "sum"

sum.abundant.UNOISE3 <- colSums(subset(reads.abundant.UNOISE3,select=-c(feature, mean.abundance)))
sum.transitional.UNOISE3 <- colSums(subset(reads.transitional.UNOISE3,select=-c(feature, mean.abundance)))
sum.rare.UNOISE3 <- colSums(subset(reads.rare.UNOISE3,select=-c(feature, mean.abundance)))
sum.UNOISE3 <- t(cbind.data.frame(sum.abundant.UNOISE3, sum.transitional.UNOISE3, sum.rare.UNOISE3))
sum.UNOISE3 <- rbind.data.frame(sum.UNOISE3, colSums(sum.UNOISE3))
row.names(sum.UNOISE3)[4] <- "sum"

sum.abundant.UPARSE <- colSums(subset(reads.abundant.UPARSE,select=-c(feature, mean.abundance)))
sum.transitional.UPARSE <- colSums(subset(reads.transitional.UPARSE,select=-c(feature, mean.abundance)))
sum.rare.UPARSE <- colSums(subset(reads.rare.UPARSE,select=-c(feature, mean.abundance)))
sum.UPARSE <- t(cbind.data.frame(sum.abundant.UPARSE, sum.transitional.UPARSE, sum.rare.UPARSE))
sum.UPARSE <- rbind.data.frame(sum.UPARSE, colSums(sum.UPARSE))
row.names(sum.UPARSE)[4] <- "sum"


# Calculate the mean and standard deviation 
AVG <- function(x)c(n=sum(!is.na(x)), mean=mean(x), sd=sd(x))

i <- m+1
j <- m+m

statistics.DaDa2.abundance <- as.data.frame(t(apply(sum.DaDa2[,c(1:m)], 1, AVG)))
statistics.DaDa2.abundance <- round(statistics.DaDa2.abundance, 4)
statistics.DaDa2.reads <- as.data.frame(t(apply(sum.DaDa2[,c(i:j)], 1, AVG)))
statistics.DaDa2.reads <- round(statistics.DaDa2.reads, 2)
AVG.abundance <- str_c(statistics.DaDa2.abundance$mean, statistics.DaDa2.abundance$sd, sep = "±")
AVG.reads <- str_c(statistics.DaDa2.reads$mean, statistics.DaDa2.reads$sd, sep = "±")
sum.DaDa2 <- cbind.data.frame(sum.DaDa2, AVG.abundance, AVG.reads)
method <- rep("DaDa2",nrow(sum.DaDa2))
sum.DaDa2 <- cbind(sum.DaDa2, method)

statistics.Deblur.abundance <- as.data.frame(t(apply(sum.Deblur[,c(1:m)], 1, AVG)))
statistics.Deblur.abundance <- round(statistics.Deblur.abundance, 4)
statistics.Deblur.reads <- as.data.frame(t(apply(sum.Deblur[,c(i:j)], 1, AVG)))
statistics.Deblur.reads <- round(statistics.Deblur.reads, 2)
AVG.abundance <- str_c(statistics.Deblur.abundance$mean, statistics.Deblur.abundance$sd, sep = "±")
AVG.reads <- str_c(statistics.Deblur.reads$mean, statistics.Deblur.reads$sd, sep = "±")
sum.Deblur <- cbind.data.frame(sum.Deblur, AVG.abundance, AVG.reads)
method <- rep("Deblur",nrow(sum.Deblur))
sum.Deblur <- cbind(sum.Deblur, method)

statistics.UNOISE3.abundance <- as.data.frame(t(apply(sum.UNOISE3[,c(1:m)], 1, AVG)))
statistics.UNOISE3.abundance <- round(statistics.UNOISE3.abundance, 4)
statistics.UNOISE3.reads <- as.data.frame(t(apply(sum.UNOISE3[,c(i:j)], 1, AVG)))
statistics.UNOISE3.reads <- round(statistics.UNOISE3.reads, 2)
AVG.abundance <- str_c(statistics.UNOISE3.abundance$mean, statistics.UNOISE3.abundance$sd, sep = "±")
AVG.reads <- str_c(statistics.UNOISE3.reads$mean, statistics.UNOISE3.reads$sd, sep = "±")
sum.UNOISE3 <- cbind.data.frame(sum.UNOISE3, AVG.abundance, AVG.reads)
method <- rep("UNOISE3",nrow(sum.UNOISE3))
sum.UNOISE3 <- cbind(sum.UNOISE3, method)

statistics.UPARSE.abundance <- as.data.frame(t(apply(sum.UPARSE[,c(1:m)], 1, AVG)))
statistics.UPARSE.abundance <- round(statistics.UPARSE.abundance, 4)
statistics.UPARSE.reads <- as.data.frame(t(apply(sum.UPARSE[,c(i:j)], 1, AVG)))
statistics.UPARSE.reads <- round(statistics.UPARSE.reads, 2)
AVG.abundance <- str_c(statistics.UPARSE.abundance$mean, statistics.UPARSE.abundance$sd, sep = "±")
AVG.reads <- str_c(statistics.UPARSE.reads$mean, statistics.UPARSE.reads$sd, sep = "±")
sum.UPARSE <- cbind.data.frame(sum.UPARSE, AVG.abundance, AVG.reads)
method <- rep("UPARSE",nrow(sum.UPARSE))
sum.UPARSE <- cbind(sum.UPARSE, method)

sum <- rbind.data.frame(sum.DaDa2, sum.Deblur, sum.UNOISE3, sum.UPARSE)

# Calculate the number of abundant, transitional and rare taxa.
n.all.DaDa2<- nrow(pick.DaDa2)
n.abundant.DaDa2 <- nrow(abundant.DaDa2)
n.transitional.DaDa2 <- nrow(transitional.DaDa2)
n.rare.DaDa2 <- nrow(rare.DaDa2)
count <- c(n.abundant.DaDa2, n.transitional.DaDa2, n.rare.DaDa2, n.all.DaDa2)
type <- c("abundant", "transitional", "rare", "sum")
count.DaDa2 <- as.data.frame(cbind(type, count))
colnames(count.DaDa2)<- c("type", "count")
method <- rep("DaDa2",nrow(count.DaDa2))
count.DaDa2 <- cbind(count.DaDa2, method )

n.all.Deblur<- nrow(pick.Deblur)
n.abundant.Deblur <- nrow(abundant.Deblur)
n.transitional.Deblur <- nrow(transitional.Deblur)
n.rare.Deblur <- nrow(rare.Deblur)
count <- c(n.abundant.Deblur, n.transitional.Deblur, n.rare.Deblur, n.all.Deblur)
type <- c("abundant", "transitional", "rare", "sum")
count.Deblur <- as.data.frame(cbind(type, count))
colnames(count.Deblur)<- c("type", "count")
method <- rep("Deblur",nrow(count.Deblur))
count.Deblur <- cbind(count.Deblur, method )

n.all.UNOISE3<- nrow(pick.UNOISE3)
n.abundant.UNOISE3 <- nrow(abundant.UNOISE3)
n.transitional.UNOISE3 <- nrow(transitional.UNOISE3)
n.rare.UNOISE3 <- nrow(rare.UNOISE3)
count <- c(n.abundant.UNOISE3, n.transitional.UNOISE3, n.rare.UNOISE3, n.all.UNOISE3)
type <- c("abundant", "transitional", "rare", "sum")
count.UNOISE3 <- as.data.frame(cbind(type, count))
colnames(count.UNOISE3)<- c("type", "count")
method <- rep("UNOISE3",nrow(count.UNOISE3))
count.UNOISE3 <- cbind(count.UNOISE3, method )

n.all.UPARSE<- nrow(pick.UPARSE)
n.abundant.UPARSE <- nrow(abundant.UPARSE)
n.transitional.UPARSE <- nrow(transitional.UPARSE)
n.rare.UPARSE <- nrow(rare.UPARSE)
count <- c(n.abundant.UPARSE, n.transitional.UPARSE, n.rare.UPARSE, n.all.UPARSE)
type <- c("abundant", "transitional", "rare", "sum")
count.UPARSE <- as.data.frame(cbind(type, count))
colnames(count.UPARSE)<- c("type", "count")
method <- rep("UPARSE",nrow(count.UPARSE))
count.UPARSE <- cbind(count.UPARSE, method )

count  <- rbind.data.frame(count.DaDa2[-4,], count.Deblur[-4,], count.UNOISE3[-4,], count.UPARSE[-4,])

# Calculate the abundance of abundant, transitional and rare taxa.
fre.abundant.DaDa2 <- sum(abundant.DaDa2$mean.abundance)
fre.transitional.DaDa2 <- sum(transitional.DaDa2$mean.abundance)
fre.rare.DaDa2 <- sum(rare.DaDa2$mean.abundance)
fre.all.DaDa2 <- fre.abundant.DaDa2+fre.transitional.DaDa2+fre.rare.DaDa2
fre.DaDa2 <- c(fre.abundant.DaDa2, fre.transitional.DaDa2, fre.rare.DaDa2, fre.all.DaDa2 )
fre.DaDa2  <- round(fre.DaDa2, 4)
fre.DaDa2 <- as.data.frame(cbind(type, fre.DaDa2))
colnames(fre.DaDa2) <- c("type", "abundance")
method <- rep("DaDa2",nrow(fre.DaDa2))
fre.DaDa2 <- cbind(fre.DaDa2, method)

fre.abundant.Deblur <- sum(abundant.Deblur$mean.abundance)
fre.transitional.Deblur <- sum(transitional.Deblur$mean.abundance)
fre.rare.Deblur <- sum(rare.Deblur$mean.abundance)
fre.all.Deblur <- fre.abundant.Deblur+fre.transitional.Deblur+fre.rare.Deblur
fre.Deblur <- c(fre.abundant.Deblur, fre.transitional.Deblur, fre.rare.Deblur, fre.all.Deblur )
fre.Deblur  <- round(fre.Deblur, 4)
fre.Deblur <- as.data.frame(cbind(type, fre.Deblur))
colnames(fre.Deblur) <- c("type", "abundance")
method <- rep("Deblur",nrow(fre.Deblur))
fre.Deblur <- cbind(fre.Deblur, method)

fre.abundant.UNOISE3 <- sum(abundant.UNOISE3$mean.abundance)
fre.transitional.UNOISE3 <- sum(transitional.UNOISE3$mean.abundance)
fre.rare.UNOISE3 <- sum(rare.UNOISE3$mean.abundance)
fre.all.UNOISE3 <- fre.abundant.UNOISE3+fre.transitional.UNOISE3+fre.rare.UNOISE3
fre.UNOISE3 <- c(fre.abundant.UNOISE3, fre.transitional.UNOISE3, fre.rare.UNOISE3, fre.all.UNOISE3 )
fre.UNOISE3  <- round(fre.UNOISE3, 4)
fre.UNOISE3 <- as.data.frame(cbind(type, fre.UNOISE3))
colnames(fre.UNOISE3) <- c("type", "abundance")
method <- rep("UNOISE3",nrow(fre.UNOISE3))
fre.UNOISE3 <- cbind(fre.UNOISE3, method)

fre.abundant.UPARSE <- sum(abundant.UPARSE$mean.abundance)
fre.transitional.UPARSE <- sum(transitional.UPARSE$mean.abundance)
fre.rare.UPARSE <- sum(rare.UPARSE$mean.abundance)
fre.all.UPARSE <- fre.abundant.UPARSE+fre.transitional.UPARSE+fre.rare.UPARSE
fre.UPARSE <- c(fre.abundant.UPARSE, fre.transitional.UPARSE, fre.rare.UPARSE, fre.all.UPARSE )
fre.UPARSE  <- round(fre.UPARSE, 4)
fre.UPARSE <- as.data.frame(cbind(type, fre.UPARSE))
colnames(fre.UPARSE) <- c("type", "abundance")
method <- rep("UPARSE",nrow(fre.UPARSE))
fre.UPARSE <- cbind(fre.UPARSE, method)

fre  <- rbind.data.frame(fre.DaDa2[-4,], fre.Deblur[-4,], fre.UNOISE3[-4,], fre.UPARSE[-4,])

# Output data.
dir.create("phylotypes")

write.table(fre, file = "phylotypes/1.Abundance.txt",col.names = NA, sep = "\t")
write.table(count, file = "phylotypes/2.Count.txt",col.names = NA, sep = "\t")
write.table(sum, file = "phylotypes/3.Details.txt",col.names = NA, sep = "\t")

write.table(abundant.DaDa2, file = "phylotypes/DaDa2_abundant_taxa.txt",col.names = NA, sep = "\t")
write.table(transitional.DaDa2, file = "phylotypes/DaDa2_transitional_taxa.txt",col.names = NA, sep = "\t")
write.table(rare.DaDa2, file = "phylotypes/DaDa2_rare_taxa.txt",col.names = NA, sep = "\t")

write.table(abundant.Deblur, file = "phylotypes/Deblur_abundant_taxa.txt",col.names = NA, sep = "\t")
write.table(transitional.Deblur, file = "phylotypes/Deblur_transitional_taxa.txt",col.names = NA, sep = "\t")
write.table(rare.Deblur, file = "phylotypes/Deblur_rare_taxa.txt",col.names = NA, sep = "\t")

write.table(abundant.UNOISE3, file = "phylotypes/UNOISE3_abundant_taxa.txt",col.names = NA, sep = "\t")
write.table(transitional.UNOISE3, file = "phylotypes/UNOISE3_transitional_taxa.txt",col.names = NA, sep = "\t")
write.table(rare.UNOISE3, file = "phylotypes/UNOISE3_rare_taxa.txt",col.names = NA, sep = "\t")

write.table(abundant.UPARSE, file = "phylotypes/UPARSE_abundant_taxa.txt",col.names = NA, sep = "\t")
write.table(transitional.UPARSE, file = "phylotypes/UPARSE_transitional_taxa.txt",col.names = NA, sep = "\t")
write.table(rare.UPARSE, file = "phylotypes/UPARSE_rare_taxa.txt",col.names = NA, sep = "\t")
