################
#Random forest 
################
library(plyr)
library(randomForest)
library(rfUtilities)

otu_table <- read.table("sequence_table.16s.filtered.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
otu_table <- t(otu_table)
metadata <- read.table("map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")

otu_nonzero_counts <- apply(otu_table, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values")

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

otu_table_rare_removed <- remove_rare(table=otu_table, cutoff_pro=0.1)
otu_table_rare_removed_norm <- sweep(otu_table_rare_removed, 2, colSums(otu_table_rare_removed) , '/')*100
otu_table_scaled <- scale(otu_table_rare_removed_norm, center = TRUE, scale = TRUE)  

#by species results
otu_table_scaled_var <- data.frame(t(otu_table_scaled))  
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Sp_by_key"]

set.seed(151)  

rf_species <- randomForest( x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)] , y=otu_table_scaled_var[ , ncol(otu_table_scaled_var)] , ntree=10000, importance=TRUE, proximities=TRUE )  

# Call:
#  randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 32.53%
# Confusion matrix:
#              gerstaeckeri indictiva lecticularia protracta sanguisuga
# gerstaeckeri           52         0            0         0          3
# indictiva               4         0            0         0          2
# lecticularia            2         0            0         0          0
# protracta               1         0            0         0          0
# sanguisuga             15         0            0         0          4
#              class.error
# gerstaeckeri  0.05454545
# indictiva     1.00000000
# lecticularia  1.00000000
# protracta     1.00000000
# sanguisuga    0.78947368

rf.significance(x=rf_species, xdata=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.012
# Model signifiant at p = 0.012
# 	 Model OOB error:  0.3253012
# 	 Random OOB error:  0.3975904
# 	 min random global error: 0.2891566
# 	 max random global error:  0.5180723
# 	 min random within class error: NA
# 	 max random within class error:  NA

#by sex results
otu_table_scaled_var <- data.frame(t(otu_table_scaled))  
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "Sex"]

set.seed(151)  

rf_sex <- randomForest( x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)] , y=otu_table_scaled_var[ , ncol(otu_table_scaled_var)] , ntree=10000, importance=TRUE, proximities=TRUE )  


# Call:
#  randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 43.37%
# Confusion matrix:
#    F  M class.error
# F 21 19   0.4750000
# M 17 26   0.3953488

rf.significance(x=rf_sex, xdata=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.159
# Model not signifiant at p = 0.159
# 	 Model OOB error:  0.4457831
# 	 Random OOB error:  0.5060241
# 	 min random global error: 0.3253012
# 	 max random global error:  0.7108434
# 	 min random within class error: 0.475
# 	 max random within class error:  0.475

#by parasite infection status results
otu_table_scaled_var <- data.frame(t(otu_table_scaled))  
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "MachRes1"]

set.seed(151)  

rf_parasite <- randomForest( x=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)] , y=otu_table_scaled_var[ , ncol(otu_table_scaled_var)] , ntree=10000, importance=TRUE, proximities=TRUE )  

# Call:
#  randomForest(x = otu_table_scaled_var[, 1:(ncol(otu_table_scaled_var) -      1)], y = otu_table_scaled_var[, ncol(otu_table_scaled_var)],      ntree = 10000, importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 37.35%
# Confusion matrix:
#          Negative Positive class.error
# Negative        9       22   0.7096774
# Positive        9       43   0.1730769

rf.significance(x=rf_parasite, xdata=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.074
# Model not signifiant at p = 0.074
# 	 Model OOB error:  0.373494
# 	 Random OOB error:  0.4457831
# 	 min random global error: 0.253012
# 	 max random global error:  0.6144578
# 	 min random within class error: 0.5806452
# 	 max random within class error:  0.5806452
