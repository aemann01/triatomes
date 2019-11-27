###########
#LIBRARIES
###########
library(tidyverse)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(compositions)
library(ape)
library(factoextra)
library(RColorBrewer)
library(philr)
library(phyloseq)
library(UpSetR)
library(reshape2)
library(vegan)
library(phylofactor)
library(ggtree)

##########
#ENV SETUP
##########
PATH="/Users/mann/github/triatomes"
setwd(PATH)
cols <- brewer.pal(6, "Set2")

###########
#LOAD DATA
###########
#load sample data
rawmetadata <- read_delim(file = file.path(PATH, "map.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries
seqtab.filtered <- read.table("sequence_table.16s.filtered.txt", header=T, row.names=1)
taxa <- read.table("assigntax/rep_set_fix_tax_assignments.txt", header=F, sep="\t", row.names=1)
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID) #record samples absent in either metadata or OTU table
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))
tree <- read.tree("rep_set.root.tre") #load representative tree

#clean up rownames in filtered sequence table - ONLY RUN THIS ONCE IF NEEDED
# new.rownames <- row.names(seqtab.filtered)
# new.rownames <- gsub("_S.*.fastq.gz", "", new.rownames)
# new.rownames <- paste("Tri", new.rownames, sep="") # append characters to numerical samples names -- these can mess things up in the long run
# rownames(seqtab.filtered) <- new.rownames
# write.table(data.frame("row_names"=rownames(seqtab.filtered),seqtab.filtered),"sequence_table.16s.filtered.txt", row.names=F, quote=F, sep="\t")
# rownames(rawmetadata) <- rawmetadata$SampleID 
# new.samplenames <- paste("Tri", rownames(rawmetadata), sep="") #add character to sample names in metadata file
# rawmetadata$SampleID <- new.samplenames

#############
#FILTER DATA
#############
#remove species with insufficient sample size
rawmetadata <- rawmetadata[rawmetadata$Sp_by_key %in% c("gerstaeckeri", "sanguisuga"),]
rownames(rawmetadata) <- rawmetadata$SampleID #set row names as the sample IDs for the metadata
seqtab.filtered <- seqtab.filtered[rownames(seqtab.filtered) %in% rownames(rawmetadata),]

################
#PHILR DISTANCE 
################
#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
ps.dat <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)
#split by species for downstream species specific analyses
ps.dat.ger <- subset_samples(ps.dat, Sp_by_key=="gerstaeckeri")
ps.dat.sang <- subset_samples(ps.dat, Sp_by_key=="sanguisuga")

#philr transform for full dataset
philr.dat <- transform_sample_counts(ps.dat, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
is.rooted(phy_tree(philr.dat)) #check that tree is rooted
is.binary.tree(phy_tree(philr.dat)) #check that multichotomies are resolved in tree
phy_tree(philr.dat) <- makeNodeLabel(phy_tree(philr.dat), method="number", prefix="n")
otu.table <- otu_table(philr.dat)
tree <- phy_tree(philr.dat)
metadata <- sample_data(philr.dat)
tax <- tax_table(philr.dat)
philr.t <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")

#philr transform gerstaeckeri only
philr.dat.ger <- transform_sample_counts(ps.dat.ger, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
phy_tree(philr.dat.ger) <- makeNodeLabel(phy_tree(philr.dat.ger), method="number", prefix="n")
otu.table <- otu_table(philr.dat.ger)
tree <- phy_tree(philr.dat.ger)
metadata <- sample_data(philr.dat.ger)
tax <- tax_table(philr.dat.ger)
philr.t.ger <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist.ger <- dist(philr.t.ger, method="euclidean")

#philr transform sanguisuga only
philr.dat.sang <- transform_sample_counts(ps.dat.sang, function(x) x+1) #add pseudocount of one to OTUs to avoid log-ratios involving zeros
phy_tree(philr.dat.sang) <- makeNodeLabel(phy_tree(philr.dat.sang), method="number", prefix="n")
otu.table <- otu_table(philr.dat.sang)
tree <- phy_tree(philr.dat.sang)
metadata <- sample_data(philr.dat.sang)
tax <- tax_table(philr.dat.sang)
philr.t.sang <- philr(otu.table, tree, part.weights="enorm.x.gm.counts", ilr.weights="blw.sqrt")
philr.dist.sang <- dist(philr.t.sang, method="euclidean")

# Heirarchical cluster dendrogram
hc <- hclust(dist(philr.t), method="complete")
df2 <- data.frame(cluster=cutree(hc,5), states=factor(hc$labels, levels=hc$labels[hc$order])) # get cluster assocaited with each sample
write.table(df2, "philr_cluster.txt", quote=F, sep="\t", col.names=NA)

hcd <- as.dendrogram(hc)
dend_data <- dendro_data(hcd, type="rectangle")
p1 <- ggplot(dend_data$segments) + geom_segment(aes(x=x,y=y, xend=xend, yend=yend)) + theme_classic() + geom_text(data = dend_data$labels, aes(x, y, label = label, hjust = 1, angle = 90)) + ylim(-2,30) + xlab("") + ylab("") + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
merge <- merge(df2, rawmetadata, by.x=c("states"), by.y=c("SampleID"))
p2 <- ggplot(merge, aes(states, y=1, fill=factor(merge$SpSex))) + geom_tile() + scale_fill_manual(values=cols) + scale_y_continuous(expand=c(0,0)) + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text=element_blank(), legend.position="none")
gp1<-ggplotGrob(p1)
gp2<-ggplotGrob(p2)
maxWidth <- grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)
pdf("figs/philr_dendrogram_spsex.pdf")
grid.arrange(gp1, gp2, ncol=1,heights=c(4/5,1/5,1/5))
dev.off()

# PCA
philr.dist <- dist(philr.t, method="euclidean")
pca <- prcomp(as.matrix(philr.dist))

pdf("figs/philr_screeplot.pdf")
screeplot(pca)
dev.off()
pdf("figs/philr_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Sp_by_key, mean.point=F, label="none", pointsize=3) + labs(title="") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-40,45) + ylim(-40,45)
dev.off()

# PCA for supplemental figure
pdf("figs/philr_bloodmeal_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Blood_meal, mean.point=F, label="none", pointsize=3) + labs(title="") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-40,45) + ylim(-40,45)
dev.off()

pdf("figs/philr_machres1_pca.pdf")
fviz_pca_ind(pca, habillage=merge$MachRes1, mean.point=F, label="none", pointsize=3) + labs(title="") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-40,45) + ylim(-40,45)
dev.off()

pdf("figs/philr_sex_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Sex, mean.point=F, label="none", pointsize=3) + labs(title="") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-40,45) + ylim(-40,45)
dev.off()

pdf("figs/philr_straintype_pca.pdf")
fviz_pca_ind(pca, habillage=merge$Strain_type, mean.point=F, label="none", pointsize=3) + labs(title="") + scale_color_brewer(palette="Dark2") + theme_minimal() + xlim(-40,45) + ylim(-40,45)
dev.off()

############
#Upset plot
############
map <- as.matrix(read.table("map.txt", header=T, sep="\t", row.names=1))
merged <- merge(seqtab.filtered, map, by="row.names")
n <- ncol(seqtab.filtered) + 1
agg <- aggregate(merged[,2:n], by=list(merge$Sp_by_key), FUN=sum)
#remove columns with only zeros
agg <- agg[,colSums(agg !=0) > 0]
rownames(agg) <- agg$Group.1
#convert to presence absence table -- ignore warnining message, still works
agg[agg>1] <- 1
#transpose again
agg <- data.frame(t(agg[,-1]))
#upsetR 
pdf("figs/upset.pdf", onefile=F)
upset(agg, order.by="freq", mainbar.y.label="Number of ASVs", sets.x.label="Shared ASVs per species", mb.ratio = c(0.55, 0.45))
dev.off()

####################
#Taxonomy barchart
####################
orderlist <- as.vector(dend_data$labels$label)
dat <- read.table("collapsed_simplified_taxonomy_for_plot.txt", header=T, sep="\t", row.names=1)
# filter out low sampled bug species
dat <- dat[,colnames(dat) %in% rownames(rawmetadata)]
# reset row names
dat <- cbind(rownames(dat), data.frame(dat, row.names=NULL))

datmelt <- melt(dat)
colnames(datmelt) <- c("taxonomy", "sampleID", "value")

datmelt$sampleID <- factor(datmelt$sampleID, levels=orderlist)
pdf("figs/taxonomy_barchart.pdf")
ggplot(datmelt, aes(fill=datmelt$taxonomy, x=datmelt$sampleID, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_classic() + theme(axis.text.x=element_text(angle=90))
dev.off()

#########
#Adonis 
#########
# is there a difference in microbial diversity across samples by some metadata category?
metadata <- as(sample_data(ps.dat), "data.frame")
metadata.ger <- as(sample_data(ps.dat.ger), "data.frame")
metadata.sang <- as(sample_data(ps.dat.sang), "data.frame")

adonis(philr.dist ~ Sex, data=metadata)
# Call:
# adonis(formula = philr.dist ~ Sex, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Sex        1     121.0 120.994  1.9118 0.02587    0.1
# Residuals 72    4556.8  63.288         0.97413
# Total     73    4677.8                 1.00000

#Ger only
adonis(philr.dist.ger ~ Sex, data=metadata.ger)

# Call:
# adonis(formula = philr.dist.ger ~ Sex, data = metadata.ger)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Sex        1      86.5  86.483  1.4273 0.02622  0.212
# Residuals 53    3211.3  60.590         0.97378
# Total     54    3297.7                 1.00000

#Sang only
adonis(philr.dist.sang ~ Sex, data=metadata.sang)

# Call:
# adonis(formula = philr.dist.sang ~ Sex, data = metadata.sang)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model     R2 Pr(>F)
# Sex        1     142.1  142.13 0.51128 0.0292  0.634
# Residuals 17    4725.7  277.98         0.9708
# Total     18    4867.9                 1.0000


adonis(philr.dist ~ Sp_by_key, data=metadata)

# Call:
# adonis(formula = philr.dist ~ Sp_by_key, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Sp_by_key  1     225.5 225.503  3.6467 0.04821  0.007 **
# Residuals 72    4452.3  61.837         0.95179
# Total     73    4677.8                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

adonis(philr.dist ~ MachRes1, data=metadata)

# Call:
# adonis(formula = philr.dist ~ MachRes1, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# MachRes1   1     185.0  184.96  2.9641 0.03954  0.036 *
# Residuals 72    4492.8   62.40         0.96046
# Total     73    4677.8                 1.00000
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#Ger only
adonis(philr.dist.ger ~ MachRes1, data=metadata.ger)

# Call:
# adonis(formula = philr.dist.ger ~ MachRes1, data = metadata.ger)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# MachRes1   1     116.2 116.202  1.9358 0.03524  0.108
# Residuals 53    3181.5  60.029         0.96476
# Total     54    3297.7                 1.00000

#Sang only
adonis(philr.dist.sang ~ MachRes1, data=metadata.sang)

# Call:
# adonis(formula = philr.dist.sang ~ MachRes1, data = metadata.sang)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# MachRes1   1     575.0  575.00   2.277 0.11812  0.108
# Residuals 17    4292.9  252.52         0.88188
# Total     18    4867.9                 1.00000

adonis(philr.dist ~ County, data=metadata)

# Call:
# adonis(formula = philr.dist ~ County, data = metadata)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# County    27    1880.6  69.650  1.1454 0.40202  0.196
# Residuals 46    2797.2  60.809         0.59798
# Total     73    4677.8                 1.00000

subset_TcI <- metadata[which(metadata$Strain_type == "TcI"),]
subset_TcIV <- metadata[which(metadata$Strain_type == "TcIV"),]
submeta <- rbind(subset_TcI, subset_TcIV)
submat <- as.matrix(philr.dist)
submat <- submat[rownames(submeta),rownames(submeta)]

adonis(submat ~ Strain_type, data=submeta)

# Call:
# adonis(formula = submat ~ Strain_type, data = submeta)

# Permutation: free
# Number of permutations: 999

# Terms added sequentially (first to last)

#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
# Strain_type  1     57.92  57.916 0.97088 0.13928   0.49
# Residuals    6    357.92  59.653         0.86072
# Total        7    415.83                 1.00000

##############
#Phylofactor
##############
OTUTable <- as.matrix(t(seqtab.filtered))
filt.list <- colnames(OTUTable)
filtmap <- rawmetadata[rawmetadata$SampleID %in% filt.list,]
filtmap <- filtmap[match(filt.list, filtmap$SampleID),]
x <- as.factor(filtmap$MachRes1) # CHANGE ME to your variable of interest
tree <- phy_tree(philr.dat)
tax <- read.table("assigntax/rep_set_tax_assignments_phylofactor.txt", sep="\t", header=T)
common.otus <- which(rowSums(OTUTable>0)>10)
OTUTable <- OTUTable[common.otus,]
tree <- ape::drop.tip(tree, setdiff(tree$tip.label, rownames(OTUTable)))
PF <- PhyloFactor(OTUTable, tree, x, nfactors=3)
PF$Data <- PF$Data[PF$tree$tip.label,]
gtree <- pf.tree(PF,layout="rectangular")
pdf("figs/phylofactor_tree.pdf")
gtree$ggplot + geom_tiplab()
dev.off()

y <- t(PF$basis[,1]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor1_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[1]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor 1') + ylim(c(-3.5,9.5))
dev.off()

wilcox.test(dat[dat$V1 == "Negative",]$V2, dat[dat$V1 == "Positive",]$V2)

# 	Wilcoxon rank sum test

# data:  dat[dat$V1 == "Negative", ]$V2 and dat[dat$V1 == "Positive", ]$V2
# W = 441, p-value = 0.02957
# alternative hypothesis: true location shift is not equal to 0

y <- t(PF$basis[,2]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor2_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[2]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor2') + ylim(c(-3.5,9.5))
dev.off()

wilcox.test(dat[dat$V1 == "Negative",]$V2, dat[dat$V1 == "Positive",]$V2)

# 	Wilcoxon rank sum test

# data:  dat[dat$V1 == "Negative", ]$V2 and dat[dat$V1 == "Positive", ]$V2
# W = 448, p-value = 0.03615
# alternative hypothesis: true location shift is not equal to 0

y <- t(PF$basis[,3]) %*% log(PF$Data)
dat <- as.data.frame(cbind(as.matrix(PF$X), (t(y))))
dat$V2 <- as.numeric(as.character(dat$V2))
pdf("figs/factor3_boxp.pdf")
ggplot(dat, aes(x=dat$V1, y=dat$V2)) + geom_boxplot(fill=gtree$legend$colors[3]) + theme_classic() + ylab("ILR abundance") + xlab("") + ggtitle('Factor3') + ylim(c(-3.5,9.5))
dev.off()

wilcox.test(dat[dat$V1 == "Negative",]$V2, dat[dat$V1 == "Positive",]$V2)

# 	Wilcoxon rank sum test with continuity correction

# data:  dat[dat$V1 == "Negative", ]$V2 and dat[dat$V1 == "Positive", ]$V2
# W = 812.5, p-value = 0.04625
# alternative hypothesis: true location shift is not equal to 0

PF$factors
#                               Group1                       Group2      ExpVar
# Factor 1 2 member Monophyletic clade 25 member Monophyletic clade 0.006195369
# Factor 2 3 member Monophyletic clade 22 member Paraphyletic clade 0.006766574
# Factor 3                         tip 21 member Paraphyletic clade 0.003837506
#                 F     Pr(>F)
# Factor 1 5.330237 0.02383237
# Factor 2 4.957810 0.02910193
# Factor 3 2.866402 0.09476980
