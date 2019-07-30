####################
#taxonomy barchart
####################
library(reshape2)
library(ggplot2)
dat <- read.table("collapsed_simplified_taxonomy.txt", header=T, sep="\t")
datmelt <- melt(dat)
pdf("taxonomy_barchart.pdf")
ggplot(datmelt, aes(fill=datmelt$SampleID, x=datmelt$variable, y=datmelt$value)) + geom_bar(stat="identity", position="fill") + theme_classic()
dev.off()


