# Calculating hamming distance.
assign.db.DaDa2.sub1$hamming <- stringdist(assign.db.DaDa2.sub1$sequences, assign.db.DaDa2.sub1$refSeqs.trim, method = "hamming")
assign.db.DaDa2.sub2$hamming <- stringdist(assign.db.DaDa2.sub2$sequences.trim, assign.db.DaDa2.sub2$refSeqs, method = "hamming")

assign.db.Deblur.sub1$hamming <- stringdist(assign.db.Deblur.sub1$sequences, assign.db.Deblur.sub1$refSeqs.trim, method = "hamming")
assign.db.Deblur.sub2$hamming <- stringdist(assign.db.Deblur.sub2$sequences.trim, assign.db.Deblur.sub2$refSeqs, method = "hamming")

assign.db.UNOISE3.sub1$hamming <- stringdist(assign.db.UNOISE3.sub1$sequences, assign.db.UNOISE3.sub1$refSeqs.trim, method = "hamming")
assign.db.UNOISE3.sub2$hamming <- stringdist(assign.db.UNOISE3.sub2$sequences.trim, assign.db.UNOISE3.sub2$refSeqs, method = "hamming")

assign.db.UPARSE.sub1$hamming <- stringdist(assign.db.UPARSE.sub1$sequences, assign.db.UPARSE.sub1$refSeqs.trim, method = "hamming")
assign.db.UPARSE.sub2$hamming <- stringdist(assign.db.UPARSE.sub2$sequences.trim, assign.db.UPARSE.sub2$refSeqs, method = "hamming")

result.DaDa2 <- rbind(subset(assign.db.DaDa2.sub1, select=-c(refSeqs.trim,length_region.trim)),subset(assign.db.DaDa2.sub2, select=-c(sequences.trim,length.trim)))
result.Deblur <- rbind(subset(assign.db.Deblur.sub1, select=-c(refSeqs.trim,length_region.trim)), subset(assign.db.Deblur.sub2, select=-c(sequences.trim,length.trim)))
result.UNOISE3 <- rbind(subset(assign.db.UNOISE3.sub1, select=-c(refSeqs.trim,length_region.trim)), subset(assign.db.UNOISE3.sub2, select=-c(sequences.trim,length.trim)))
result.UPARSE <- rbind(subset(assign.db.UPARSE.sub1, select=-c(refSeqs.trim,length_region.trim)), subset(assign.db.UPARSE.sub2, select=-c(sequences.trim,length.trim)))

DaDa2 <- split(result.DaDa2, result.DaDa2$feature)
hamming.DaDa2 <- as.data.frame(sapply(DaDa2,FUN=function(x) min(x$hamming)))
colnames(hamming.DaDa2) <- "hamming"
feature.DaDa2 <-row.names(hamming.DaDa2)
hamming.DaDa2 <- cbind(feature.DaDa2, hamming.DaDa2)

Deblur <- split(result.Deblur, result.Deblur$feature)
hamming.Deblur <- as.data.frame(sapply(Deblur,FUN=function(x) min(x$hamming)))
colnames(hamming.Deblur) <- "hamming"
feature.Deblur <-row.names(hamming.Deblur)
hamming.Deblur <- cbind(feature.Deblur, hamming.Deblur)

UNOISE3 <- split(result.UNOISE3, result.UNOISE3$feature)
hamming.UNOISE3 <- as.data.frame(sapply(UNOISE3,FUN=function(x) min(x$hamming)))
colnames(hamming.UNOISE3) <- "hamming"
feature.UNOISE3 <-row.names(hamming.UNOISE3)
hamming.UNOISE3 <- cbind(feature.UNOISE3, hamming.UNOISE3)

UPARSE <- split(result.UPARSE, result.UPARSE$feature)
hamming.UPARSE <- as.data.frame(sapply(UPARSE,FUN=function(x) min(x$hamming)))
colnames(hamming.UPARSE) <- "hamming"
feature.UPARSE <-row.names(hamming.UPARSE)
hamming.UPARSE <- cbind(feature.UPARSE, hamming.UPARSE)

hamming.DaDa2 <-  merge(seq.DaDa2, hamming.DaDa2, by.x = "feature.DaDa2", by.y = "feature.DaDa2") 
hamming.Deblur <-  merge(seq.Deblur, hamming.Deblur, by.x = "feature.Deblur", by.y = "feature.Deblur") 
hamming.UNOISE3 <-  merge(seq.UNOISE3, hamming.UNOISE3, by.x = "feature.UNOISE3", by.y = "feature.UNOISE3") 
hamming.UPARSE <-  merge(seq.UPARSE, hamming.UPARSE, by.x = "feature.UPARSE", by.y = "feature.UPARSE") 

# Statistics.
hamming.DaDa2.all <- nrow(hamming.DaDa2)
hamming.DaDa2.HD0 <- nrow(hamming.DaDa2[which(hamming.DaDa2$hamming == 0),])
hamming.DaDa2.HD0.ratio <- round(hamming.DaDa2.HD0/nrow(hamming.DaDa2)*100, 2)
hamming.DaDa2.HD1 <- nrow(hamming.DaDa2[which(hamming.DaDa2$hamming == 1 ),])
hamming.DaDa2.HD1.ratio <- round(hamming.DaDa2.HD1/hamming.DaDa2.all*100, 2)
hamming.DaDa2.HD2_ <- nrow(hamming.DaDa2[which(hamming.DaDa2$hamming>= 2),])
hamming.DaDa2.HD2_.ratio <- round(hamming.DaDa2.HD2_/hamming.DaDa2.all*100, 2)
hamming.DaDa2.summary <- round(summary(hamming.DaDa2$hamming), 2)
algorithm <- "DaDa2"
hamming.DaDa2.abstract <- cbind(hamming.DaDa2.all, hamming.DaDa2.HD0.ratio,
                               hamming.DaDa2.HD1.ratio, hamming.DaDa2.HD2_.ratio,
                               t(hamming.DaDa2.summary), algorithm)
colnames(hamming.DaDa2.abstract) <- c("features", "HD0", "HD1", ">=HD2", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max","algorithm")

hamming.Deblur.all <- nrow(hamming.Deblur)
hamming.Deblur.HD0 <- nrow(hamming.Deblur[which(hamming.Deblur$hamming == 0),])
hamming.Deblur.HD0.ratio <- round(hamming.Deblur.HD0/nrow(hamming.Deblur)*100, 2)
hamming.Deblur.HD1 <- nrow(hamming.Deblur[which(hamming.Deblur$hamming == 1 ),])
hamming.Deblur.HD1.ratio <- round(hamming.Deblur.HD1/hamming.Deblur.all*100, 2)
hamming.Deblur.HD2_ <- nrow(hamming.Deblur[which(hamming.Deblur$hamming>= 2),])
hamming.Deblur.HD2_.ratio <- round(hamming.Deblur.HD2_/hamming.Deblur.all*100, 2)
hamming.Deblur.summary <- round(summary(hamming.Deblur$hamming), 2)
algorithm <- "Deblur"
hamming.Deblur.abstract <- cbind(hamming.Deblur.all, hamming.Deblur.HD0.ratio,
                                hamming.Deblur.HD1.ratio, hamming.Deblur.HD2_.ratio,
                                t(hamming.Deblur.summary), algorithm)
colnames(hamming.Deblur.abstract) <- c("features", "HD0", "HD1", ">=HD2", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max","algorithm")

hamming.UNOISE3.all <- nrow(hamming.UNOISE3)
hamming.UNOISE3.HD0 <- nrow(hamming.UNOISE3[which(hamming.UNOISE3$hamming == 0),])
hamming.UNOISE3.HD0.ratio <- round(hamming.UNOISE3.HD0/nrow(hamming.UNOISE3)*100, 2)
hamming.UNOISE3.HD1 <- nrow(hamming.UNOISE3[which(hamming.UNOISE3$hamming == 1 ),])
hamming.UNOISE3.HD1.ratio <- round(hamming.UNOISE3.HD1/hamming.UNOISE3.all*100, 2)
hamming.UNOISE3.HD2_ <- nrow(hamming.UNOISE3[which(hamming.UNOISE3$hamming>= 2),])
hamming.UNOISE3.HD2_.ratio <- round(hamming.UNOISE3.HD2_/hamming.UNOISE3.all*100, 2)
hamming.UNOISE3.summary <- round(summary(hamming.UNOISE3$hamming), 2)
algorithm <- "UNOISE3"
hamming.UNOISE3.abstract <- cbind(hamming.UNOISE3.all, hamming.UNOISE3.HD0.ratio,
                                hamming.UNOISE3.HD1.ratio, hamming.UNOISE3.HD2_.ratio,
                                t(hamming.UNOISE3.summary), algorithm)
colnames(hamming.UNOISE3.abstract) <- c("features", "HD0", "HD1", ">=HD2", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max","algorithm")

hamming.UPARSE.all <- nrow(hamming.UPARSE)
hamming.UPARSE.HD0 <- nrow(hamming.UPARSE[which(hamming.UPARSE$hamming == 0),])
hamming.UPARSE.HD0.ratio <- round(hamming.UPARSE.HD0/nrow(hamming.UPARSE)*100, 2)
hamming.UPARSE.HD1 <- nrow(hamming.UPARSE[which(hamming.UPARSE$hamming == 1 ),])
hamming.UPARSE.HD1.ratio <- round(hamming.UPARSE.HD1/hamming.UPARSE.all*100, 2)
hamming.UPARSE.HD2_ <- nrow(hamming.UPARSE[which(hamming.UPARSE$hamming>= 2),])
hamming.UPARSE.HD2_.ratio <- round(hamming.UPARSE.HD2_/hamming.UPARSE.all*100, 2)
hamming.UPARSE.summary <- round(summary(hamming.UPARSE$hamming), 2)
algorithm <- "UPARSE"
hamming.UPARSE.abstract <- cbind(hamming.UPARSE.all,  hamming.UPARSE.HD0.ratio,
                                hamming.UPARSE.HD1.ratio, hamming.UPARSE.HD2_.ratio,
                                t(hamming.UPARSE.summary), algorithm)
colnames(hamming.UPARSE.abstract) <- c("features", "HD0", "HD1", ">=HD2", "Min", "1st_Qu", "Median", "Mean", "3rd_Qu", "Max","algorithm")

result.abstract <- as.data.frame(rbind(hamming.DaDa2.abstract, hamming.Deblur.abstract, hamming.UNOISE3.abstract, hamming.UPARSE.abstract))
rownames(result.abstract) <- result.abstract$algorithm

# Illustrate figure.
plot.abstract <- result.abstract[,c(2:4,11)]
colnames(plot.abstract) <- c("0", "1", ">=2", "algorithm")
plot.abstract <- melt(plot.abstract, id.vars = "algorithm", variable.name="hamming")

colors <- c("#C2262C", "#557EB5", "#75AD51", "#BBBBBB")
scaleFUN <- function(x){sprintf("%.1f", x)}
P <- ggplot(plot.abstract, aes(x = hamming, y = as.numeric(value), group = algorithm)) +
  geom_point(aes(colour = algorithm), size = 3, alpha = 0.5) +
  scale_color_manual(values = colors) +
  scale_y_continuous(labels = scaleFUN) +
  facet_wrap(~hamming, scales="free") +
  theme_bw() +
  labs(x = "Hamming", y = "Percentage(%)") +
  theme(text = element_text(size =10),
        legend.position="top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        strip.text.x = element_blank()
        ) 
P

# Output data.
dir.create("hamming")

write.table(hamming.DaDa2, file ="hamming/hamming.DaDa2.txt")
write.table(hamming.Deblur, file ="hamming/hamming.Deblur.txt")
write.table(hamming.UNOISE3, file ="hamming/hamming.UNOISE3.txt")
write.table(hamming.UPARSE, file ="hamming/hamming.UPARSE.txt")
write.table(result.abstract, file = "hamming/Hamming_distance_result.txt",col.names = NA, sep = "\t")
ggsave("hamming/Hamming_distance_percentage.pdf", P, width = 150, height = 120, units = 'mm' )

