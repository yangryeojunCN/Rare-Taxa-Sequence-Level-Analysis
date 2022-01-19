# Delete taxonomy classification  in Species level.
taxa_region <-row.names(ref_region)
ref_region <- cbind(taxa_region, ref_region)
colnames(ref_region) <- c("taxa", "refSeqs")
length_region <- nchar(ref_region$refSeqs)
ref_region <- cbind(ref_region, length_region)

region <-merge(ref_region, tax_region, by.x = "taxa", by.y = "Feature.ID") 
region.sep <- separate(data = region, col = Taxon, 
                   into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")
region.sep <- region.sep[, -10]

region.Taxon  <- str_c(region.sep$Kingdom, region.sep$Phylum, region.sep$Class, region.sep$Order, region.sep$Family,  region.sep$Genus, sep = ";")
region.sub <- cbind(region.sep[,c(1,2,3)], region.Taxon)
colnames(region.sub)[4] <- "Taxon"

# Match the representative sequences of rare taxa. 
seq.rare.DaDa2 <-merge(rare.DaDa2, seq.DaDa2, by.x = "feature", by.y = "feature.DaDa2") 
assign.DaDa2 <-merge(seq.rare.DaDa2, tax.DaDa2, by.x = "feature", by.y = "Feature.ID") 

seq.rare.Deblur <-merge(rare.Deblur, seq.Deblur, by.x = "feature", by.y = "feature.Deblur") 
assign.Deblur <-merge(seq.rare.Deblur, tax.Deblur, by.x = "feature", by.y = "Feature.ID") 

seq.rare.UNOISE3 <-merge(rare.UNOISE3, seq.UNOISE3, by.x = "feature", by.y = "feature.UNOISE3") 
assign.UNOISE3 <-merge(seq.rare.UNOISE3, tax.UNOISE3, by.x = "feature", by.y = "Feature.ID") 

seq.rare.UPARSE <-merge(rare.UPARSE, seq.UPARSE, by.x = "feature", by.y = "feature.UPARSE") 
assign.UPARSE <-merge(seq.rare.UPARSE, tax.UPARSE, by.x = "feature", by.y = "Feature.ID") 

# Match the sequences from database.
assign.db.DaDa2 <- merge(assign.DaDa2, region, by = "Taxon", all.x = TRUE) 
assign.db.DaDa2.filter <- filter(assign.db.DaDa2, !is.na(assign.db.DaDa2$refSeqs))
assign.db.DaDa2[is.na(assign.db.DaDa2)]<- "unassigned"
assign.db.DaDa2.Na <- filter(assign.db.DaDa2, grepl("unassigned", taxa))
assign.db.DaDa2.Na <- filter(assign.db.DaDa2.Na, grepl("g_", Taxon))
assign.db.DaDa2.Na <-subset(assign.db.DaDa2.Na, select=-c(taxa,refSeqs,length_region))
assign.db.DaDa2.Na <-merge(assign.db.DaDa2.Na, region.sub, by = "Taxon", all.x = TRUE)
assign.db.DaDa2 <-rbind(assign.db.DaDa2.filter, assign.db.DaDa2.Na)

assign.db.Deblur <- merge(assign.Deblur, region, by = "Taxon", all.x = TRUE) 
assign.db.Deblur.filter <- filter(assign.db.Deblur, !is.na(assign.db.Deblur$refSeqs))
assign.db.Deblur[is.na(assign.db.Deblur)]<- "unassigned"
assign.db.Deblur.Na <- filter(assign.db.Deblur, grepl("unassigned", taxa))
assign.db.Deblur.Na <- filter(assign.db.Deblur.Na, grepl("g_", Taxon))
assign.db.Deblur.Na <- subset(assign.db.Deblur.Na, select=-c(taxa,refSeqs,length_region))
assign.db.Deblur.Na <- merge(assign.db.Deblur.Na, region.sub, by = "Taxon", all.x = TRUE)
assign.db.Deblur <- rbind(assign.db.Deblur.filter, assign.db.Deblur.Na)

assign.db.UNOISE3 <- merge(assign.UNOISE3, region, by = "Taxon", all.x = TRUE) 
assign.db.UNOISE3.filter <- filter(assign.db.UNOISE3, !is.na(assign.db.UNOISE3$refSeqs))
assign.db.UNOISE3[is.na(assign.db.UNOISE3)]<- "unassigned"
assign.db.UNOISE3.Na <- filter(assign.db.UNOISE3, grepl("unassigned", taxa))
assign.db.UNOISE3.Na <- filter(assign.db.UNOISE3.Na, grepl("g_", Taxon))
assign.db.UNOISE3.Na <- subset(assign.db.UNOISE3.Na, select=-c(taxa,refSeqs,length_region))
assign.db.UNOISE3.Na <- merge(assign.db.UNOISE3.Na, region.sub, by = "Taxon", all.x = TRUE)
assign.db.UNOISE3 <- rbind(assign.db.UNOISE3.filter, assign.db.UNOISE3.Na)

assign.db.UPARSE <- merge(assign.UPARSE, region, by = "Taxon", all.x = TRUE) 
assign.db.UPARSE.filter <- filter(assign.db.UPARSE, !is.na(assign.db.UPARSE$refSeqs))
assign.db.UPARSE[is.na(assign.db.UPARSE)]<- "unassigned"
assign.db.UPARSE.Na <- filter(assign.db.UPARSE, grepl("unassigned", taxa))
assign.db.UPARSE.Na <- filter(assign.db.UPARSE.Na, grepl("g_", Taxon))
assign.db.UPARSE.Na <- subset(assign.db.UPARSE.Na, select=-c(taxa,refSeqs,length_region))
assign.db.UPARSE.Na <-merge(assign.db.UPARSE.Na, region.sub, by = "Taxon", all.x = TRUE)
assign.db.UPARSE <-rbind(assign.db.UPARSE.filter, assign.db.UPARSE.Na)

# Trim the length of sequences from database.
assign.db.DaDa2$diff <- assign.db.DaDa2$length_region - assign.db.DaDa2$length
assign.db.Deblur$diff <- assign.db.Deblur$length_region - assign.db.Deblur$length
assign.db.UNOISE3$diff <- assign.db.UNOISE3$length_region - assign.db.UNOISE3$length
assign.db.UPARSE$diff <- assign.db.UPARSE$length_region - assign.db.UPARSE$length

assign.db.DaDa2.sub1 <- assign.db.DaDa2[which(assign.db.DaDa2$diff >= 0),]
assign.db.DaDa2.sub2 <- assign.db.DaDa2[which(assign.db.DaDa2$diff < 0),]
assign.db.DaDa2.sub1$refSeqs.trim <- substring(assign.db.DaDa2.sub1$refSeqs, 1, assign.db.DaDa2.sub1$length)
assign.db.DaDa2.sub1$length_region.trim <- nchar(assign.db.DaDa2.sub1$refSeqs.trim)
assign.db.DaDa2.sub2$sequences.trim <- substring(assign.db.DaDa2.sub2$sequence, 1, assign.db.DaDa2.sub2$length_region)
assign.db.DaDa2.sub2$length.trim <- nchar(assign.db.DaDa2.sub2$sequences.trim)

assign.db.Deblur.sub1 <- assign.db.Deblur[which(assign.db.Deblur$diff >= 0),]
assign.db.Deblur.sub2 <- assign.db.Deblur[which(assign.db.Deblur$diff < 0),]
assign.db.Deblur.sub1$refSeqs.trim <- substring(assign.db.Deblur.sub1$refSeqs, 1, assign.db.Deblur.sub1$length)
assign.db.Deblur.sub1$length_region.trim <- nchar(assign.db.Deblur.sub1$refSeqs.trim)
assign.db.Deblur.sub2$sequences.trim <- substring(assign.db.Deblur.sub2$sequence, 1, assign.db.Deblur.sub2$length_region)
assign.db.Deblur.sub2$length.trim <- nchar(assign.db.Deblur.sub2$sequences.trim)

assign.db.UNOISE3.sub1 <- assign.db.UNOISE3[which(assign.db.UNOISE3$diff >= 0),]
assign.db.UNOISE3.sub2 <- assign.db.UNOISE3[which(assign.db.UNOISE3$diff < 0),]
assign.db.UNOISE3.sub1$refSeqs.trim <- substring(assign.db.UNOISE3.sub1$refSeqs, 1, assign.db.UNOISE3.sub1$length)
assign.db.UNOISE3.sub1$length_region.trim <- nchar(assign.db.UNOISE3.sub1$refSeqs.trim)
assign.db.UNOISE3.sub2$sequences.trim <- substring(assign.db.UNOISE3.sub2$sequence, 1, assign.db.UNOISE3.sub2$length_region)
assign.db.UNOISE3.sub2$length.trim <- nchar(assign.db.UNOISE3.sub2$sequences.trim)

assign.db.UPARSE.sub1 <- assign.db.UPARSE[which(assign.db.UPARSE$diff >= 0),]
assign.db.UPARSE.sub2 <- assign.db.UPARSE[which(assign.db.UPARSE$diff < 0),]
assign.db.UPARSE.sub1$refSeqs.trim <- substring(assign.db.UPARSE.sub1$refSeqs, 1, assign.db.UPARSE.sub1$length)
assign.db.UPARSE.sub1$length_region.trim <- nchar(assign.db.UPARSE.sub1$refSeqs.trim)
assign.db.UPARSE.sub2$sequences.trim <- substring(assign.db.UPARSE.sub2$sequence, 1, assign.db.UPARSE.sub2$length_region)
assign.db.UPARSE.sub2$length.trim <- nchar(assign.db.UPARSE.sub2$sequences.trim)

# Output data.
dir.create("assigned")

write.table(assign.db.DaDa2.sub1, file ="assigned/assign.db.DaDa2.sub1.txt")
write.table(assign.db.DaDa2.sub2, file ="assigned/assign.db.DaDa2.sub2.txt")
write.table(assign.db.Deblur.sub1, file ="assigned/assign.db.Deblur.sub1.txt")
write.table(assign.db.Deblur.sub2, file ="assigned/assign.db.Deblur.sub2.txt")
write.table(assign.db.UNOISE3.sub1, file ="assigned/assign.db.UNOISE3.sub1.txt")
write.table(assign.db.UNOISE3.sub2, file ="assigned/assign.db.UNOISE3.sub2.txt")
write.table(assign.db.UPARSE.sub1, file ="assigned/assign.db.UPARSE.sub1.txt")
write.table(assign.db.UPARSE.sub2, file ="assigned/assign.db.UPARSE.sub2.txt")