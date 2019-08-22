ger <- rawmetadata[rawmetadata$Sp_by_key == "gerstaeckeri",]$SampleID
ger.subset <- subset(seqtab.filtered, rownames(seqtab.filtered) %in% ger)

metadata.sub <- subset(rawmetadata, rawmetadata$SampleID %in% ger)

rownames(metadata.sub) <- metadata.sub$SampleID

ps.dat <- phyloseq(otu_table(ger.subset, taxa_are_rows=FALSE), 
                          sample_data(metadata.sub), 
                          tax_table(as.matrix(taxa[1])), tree)

philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")

philr.dist <- dist(philr.t, method="euclidean")

# is gerstaeckeri significant by itself?
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ MachRes1, data=metadata)

# Call:
# adonis(formula = philr.dist ~ MachRes1, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# MachRes1   1     116.2 116.202  1.9358 0.03524  0.128
# Residuals 53    3181.5  60.029         0.96476
# Total     54    3297.7                 1.00000

sang <- rawmetadata[rawmetadata$Sp_by_key == "sanguisuga",]$SampleID
sang.subset <- subset(seqtab.filtered, rownames(seqtab.filtered) %in% sang)

metadata.sub <- subset(rawmetadata, rawmetadata$SampleID %in% sang.subset)

rownames(metadata.sub) <- metadata.sub$SampleID

ps.dat <- phyloseq(otu_table(sang.subset, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)

philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")

philr.dist <- dist(philr.t, method="euclidean")

# is sanguisuga significant by itself?
metadata <- as(sample_data(ps.dat), "data.frame")
adonis(philr.dist ~ MachRes1, data=metadata)

# Call:
# adonis(formula = philr.dist ~ MachRes1, data = metadata)

# Permutation: free

# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# MachRes1   1     575.0  575.00   2.277 0.11812  0.109
# Residuals 17    4292.9  252.52         0.88188
# Total     18    4867.9                 1.00000