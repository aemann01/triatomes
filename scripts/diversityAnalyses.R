####Diversity analyses (Weighted Unifrac, Aitchison Distance, Bray Curtis Dissimilarity)

###########
#LIBRARIES
###########
library(phylofactor)
library(ape)
library(phyloseq)
library(tidyverse)
library(reshape2)
library(stringr)
library(data.table)
library(broom)
library(qualpalr)
library(viridis)
library(ShortRead)
library(Biostrings)
library(seqinr)
library(compositions)
library(ape)
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(dendextend)
library(phangorn)
library(vegan)
library(factoextra)
library(ranacapa)

PATH="/Users/mann/github/triatomes"

setwd(PATH)

#########################
#LOAD, FILTER, PREP DATA
#########################
#load sample data
rawmetadata <- read_delim(file = file.path(PATH, "map.txt"), # file.path() is used for cross-platform compatibility
                          "\t", # the text file is tab delimited
                          escape_double = FALSE, # the imported text file does not 'escape' quotation marks by wrapping them with more quotation marks
                          trim_ws = TRUE) # remove leading and trailing spaces from character string entries

seqtab.filtered <- read.table("sequence_table.16s.filtered.gsOnly.txt", header=T, row.names=1)
taxa <- read.table("assigntax/rep_set_fix_tax_assignments.txt", header=F, sep="\t", row.names=1)

#record samples absent in either metadata or OTU table
notinmeta <- setdiff(row.names(seqtab.filtered), rawmetadata$SampleID)
notinraw <- setdiff(rawmetadata$SampleID, row.names(seqtab.filtered))

#load representative tree
tree <- read.tree("rep_set.filt.tre")

#clean up rownames in filtered sequence table
# new.rownames <- row.names(seqtab.filtered)
# new.rownames <- gsub("_S.*.fastq.gz", "", new.rownames)
# new.rownames <- paste("Tri", new.rownames, sep="") # append characters to numerical samples names -- these can mess things up in the long run
# rownames(seqtab.filtered) <- new.rownames
# write.table(data.frame("row_names"=rownames(seqtab.filtered),seqtab.filtered),"sequence_table.16s.filtered.txt", row.names=F, quote=F, sep="\t")

#add character to sample names in metadata file
#set row names as the sample IDs for the metadata
# rownames(rawmetadata) <- rawmetadata$SampleID
# new.samplenames <- paste("Tri", rownames(rawmetadata), sep="")
# rawmetadata$SampleID <- new.samplenames

rownames(rawmetadata) <- rawmetadata$SampleID
#create phyloseq object from "seqtab.filtered", "rawmetadata", "taxa", and "tree"
ps.dada2_join <- phyloseq(otu_table(seqtab.filtered, taxa_are_rows=FALSE), 
                          sample_data(rawmetadata), 
                          tax_table(as.matrix(taxa[1])), tree)

###########################
#BRAY CURTIS DISSIMILARITY
###########################
ps.rare <- rarefy_even_depth(ps.dada2_join, rngseed=2867)
bray_pcoa <- ordinate(physeq=ps.rare, method="PCoA", distance="bray")
pdf("figs/braycurtis_pca_update.pdf")
plot_ordination(ps.rare, bray_pcoa, color="Sp_by_key", shape="MachRes1") + theme_minimal()
dev.off()

##########################
#WEIGHTED UNIFRAC DISTANCE
##########################
wuni_pcoa <- ordinate(physeq=ps.rare, method="PCoA", distance="wunifrac")
pdf("figs/wunifrac_pca_update.pdf")
plot_ordination(ps.rare, wuni_pcoa, color="Sp_by_key", shape="MachRes1") + theme_minimal()
dev.off()

########################
#RAREFACTION CURVE PLOT
########################
pdf("figs/rarefaction_curve.pdf")
ggrare(ps.dada2_join, color="Sp_by_key") + theme_minimal()
dev.off()

########################
#ALPHA DIVERSITY PLOT
########################
pdf("figs/adiv_spbykey.pdf")
plot_richness(ps.rare, measures=c("Observed", "Shannon"), x="Sp_by_key") + theme_minimal()
dev.off()

pdf("figs/adiv_machres1.pdf")
plot_richness(ps.rare, measures=c("Observed", "Shannon"), x="MachRes1") + theme_minimal()
dev.off()

pdf("figs/adiv.pdf")
plot_richness(ps.rare, measures=c("Observed", "Shannon"), x="Sp_by_key", color="MachRes1") + theme_minimal()
dev.off()

# signficance test
p <- estimate_richness(ps.rare, measures=c("Observed", "Shannon"))
d <- sample_data(ps.rare)
sp_ger <- p[d[,"Sp_by_key"] == "gerstaeckeri",]
sp_san <- p[d[,"Sp_by_key"] == "sanguisuga",]

wilcox.test(sp_ger$Observed, sp_san$Observed)

# 	Wilcoxon rank sum test with continuity correction

# data:  sp_ger$Observed and sp_san$Observed
# W = 346.5, p-value = 0.02959
# alternative hypothesis: true location shift is not equal to 0

#remove T. sanguisuga outlier -- significant still?
wilcox.test(sp_ger$Observed, sp_san$Observed)

# 	Wilcoxon rank sum test with continuity correction

# data:  sp_ger$Observed and sp_san$Observed
# W = 346.5, p-value = 0.05772
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(sp_ger$Shannon, sp_san$Shannon)

# 	Wilcoxon rank sum test with continuity correction

# data:  sp_ger$Shannon and sp_san$Shannon
# W = 502, p-value = 0.8045
# alternative hypothesis: true location shift is not equal to 0

cruz_pos <- p[d[,"MachRes1"] == "Positive",]
cruz_neg <- p[d[,"MachRes1"] == "Negative",]

wilcox.test(cruz_pos$Observed, cruz_neg$Observed)

# 	Wilcoxon rank sum test with continuity correction

# data:  cruz_pos$Observed and cruz_neg$Observed
# W = 670.5, p-value = 0.6896
# alternative hypothesis: true location shift is not equal to 0

wilcox.test(cruz_pos$Shannon, cruz_neg$Shannon)

# 	Wilcoxon rank sum test

# data:  cruz_pos$Shannon and cruz_neg$Shannon
# W = 647, p-value = 0.8935
# alternative hypothesis: true location shift is not equal to 0





