################
#Random forest 
################
library(plyr)
library(randomForest)
library(rfUtilities)

otu_table <- read.table("sequence_table.16s.filtered.gsOnly.txt", sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
otu_table <- t(otu_table)
metadata <- read.table("map.txt", sep="\t", header=T, row.names=1, stringsAsFactors=TRUE, comment.char="")
metadata <- metadata[metadata$Sp_by_key %in% c("gerstaeckeri", "sanguisuga"),]
metadata$Sp_by_key <- factor(metadata$Sp_by_key)

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

#         OOB estimate of  error rate: 24.32%
# Confusion matrix:
#              gerstaeckeri sanguisuga class.error
# gerstaeckeri           51          4  0.07272727
# sanguisuga             14          5  0.73684211

rf.significance(x=rf_species, xdata=otu_table_scaled_var[,1:(ncol(otu_table_scaled_var)-1)], nperm=1000, ntree=1000)

# Number of permutations:  1000
# p-value:  0.022
# Model signifiant at p = 0.022
#    Model OOB error:  0.2432432
#    Random OOB error:  0.3108108
#    min random global error: 0.2027027
#    max random global error:  0.4324324
#    min random within class error: 0.6315789
#    max random within class error:  0.6315789

#even sampling depth
otu_table_ger <- subset(otu_table_scaled_var, var=="gerstaeckeri")
otu_table_san <- subset(otu_table_scaled_var, var=="sanguisuga")
#randomly sample 19 rows
randomRows = function(df,n){
   return(df[sample(nrow(df),n),])
}
otu_table_ger_sub <- randomRows(otu_table_ger, 19)
otu_table_san_sub <- randomRows(otu_table_san, 19)
#bind together
subset_otu_table <- rbind(otu_table_ger_sub, otu_table_san_sub)
#reset factors
subset_otu_table$var <- factor(subset_otu_table$var)
# run random forest
set.seed(151)  
rf_species <- randomForest( x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  

# Call:
#  randomForest(x = subset_otu_table[, 1:(ncol(subset_otu_table) -      1)], y = subset_otu_table[, ncol(subset_otu_table)], ntree = 10000,      importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 34.21%
# Confusion matrix:
#              gerstaeckeri sanguisuga class.error
# gerstaeckeri           14          5   0.2631579
# sanguisuga              8         11   0.4210526

rf.significance(x=rf_species, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.034
# Model signifiant at p = 0.034
#    Model OOB error:  0.3421053
#    Random OOB error:  0.5263158
#    min random global error: 0.2105263
#    max random global error:  0.8421053
#    min random within class error: 0.4210526
#    max random within class error:  0.4210526

# T. gerstaeckeri by sex
otu_table_scaled_var <- data.frame(t(otu_table_scaled))  
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "SpSex"]
otu_table_gerM <- subset(otu_table_scaled_var, var=="gerstaeckeri_M")
otu_table_gerF <- subset(otu_table_scaled_var, var=="gerstaeckeri_F")

#randomly select 26 from each
otu_table_gerM_sub <- randomRows(otu_table_gerM, 26)
otu_table_gerF_sub <- randomRows(otu_table_gerF, 26)

#bind together
subset_otu_table <- rbind(otu_table_gerM_sub, otu_table_gerF_sub)
#reset factors
subset_otu_table$var <- factor(subset_otu_table$var)
# run random forest
set.seed(151) 
rf_sex <- randomForest( x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  

# Call:
#  randomForest(x = subset_otu_table[, 1:(ncol(subset_otu_table) -      1)], y = subset_otu_table[, ncol(subset_otu_table)], ntree = 10000,      importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 53.85%
# Confusion matrix:
#                gerstaeckeri_F gerstaeckeri_M class.error
# gerstaeckeri_F             13             13   0.5000000
# gerstaeckeri_M             15             11   0.5769231


rf.significance(x=rf_sex, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.566
# Model not signifiant at p = 0.566
#    Model OOB error:  0.5384615
#    Random OOB error:  0.5192308
#    min random global error: 0.25
#    max random global error:  0.7884615
#    min random within class error: 0.4444444
#    max random within class error:  0.4444444

# T. sanguisuga by sex
otu_table_gerM <- subset(otu_table_scaled_var, var=="sanguisuga_M")
otu_table_gerF <- subset(otu_table_scaled_var, var=="sanguisuga_F")

#randomly select 9 from each
otu_table_gerM_sub <- randomRows(otu_table_gerM, 9)
otu_table_gerF_sub <- randomRows(otu_table_gerF, 9)

#bind together
subset_otu_table <- rbind(otu_table_gerM_sub, otu_table_gerF_sub)
#reset factors
subset_otu_table$var <- factor(subset_otu_table$var)
# run random forest
set.seed(151) 
rf_sex <- randomForest( x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  

# Call:
#  randomForest(x = subset_otu_table[, 1:(ncol(subset_otu_table) -      1)], y = subset_otu_table[, ncol(subset_otu_table)], ntree = 10000,      importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 33.33%
# Confusion matrix:
#              sanguisuga_F sanguisuga_M class.error
# sanguisuga_F            7            2   0.2222222
# sanguisuga_M            4            5   0.4444444

rf.significance(x=rf_sex, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.053
# Model not signifiant at p = 0.053
#    Model OOB error:  0.3333333
#    Random OOB error:  0.5555556
#    min random global error: 0.1111111
#    max random global error:  1
#    min random within class error: NA
#    max random within class error:  NA

# T. gerstaeckeri by parasite infection status
metadata$sp_machres <- paste(metadata$Sp_by_key, metadata$MachRes1, sep="")
otu_table_scaled_var <- data.frame(t(otu_table_scaled)) 
otu_table_scaled_var$var <- metadata[rownames(otu_table_scaled_var), "sp_machres"]
otu_table_gerM <- subset(otu_table_scaled_var, var=="gerstaeckeriPositive")
otu_table_gerF <- subset(otu_table_scaled_var, var=="gerstaeckeriNegative")

#randomly select 18 from each
otu_table_gerM_sub <- randomRows(otu_table_gerM, 18)
otu_table_gerF_sub <- randomRows(otu_table_gerF, 18)

#bind together
subset_otu_table <- rbind(otu_table_gerM_sub, otu_table_gerF_sub)
#reset factors
subset_otu_table$var <- factor(subset_otu_table$var)
# run random forest
set.seed(151) 
rf_par <- randomForest( x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  

# Call:
#  randomForest(x = subset_otu_table[, 1:(ncol(subset_otu_table) -      1)], y = subset_otu_table[, ncol(subset_otu_table)], ntree = 10000,      importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 47.22%
# Confusion matrix:
#                      gerstaeckeriNegative gerstaeckeriPositive class.error
# gerstaeckeriNegative                    8                   10   0.5555556
# gerstaeckeriPositive                    7                   11   0.3888889


rf.significance(x=rf_par, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.258
# Model not signifiant at p = 0.258
#    Model OOB error:  0.4722222
#    Random OOB error:  0.5277778
#    min random global error: 0.1944444
#    max random global error:  0.8333333
#    min random within class error: 0.4666667
#    max random within class error:  0.4666667

# T. sanguisuga by parasite status
otu_table_gerM <- subset(otu_table_scaled_var, var=="sanguisugaPositive")
otu_table_gerF <- subset(otu_table_scaled_var, var=="sanguisugaNegative")

#randomly select 9 from each
otu_table_gerM_sub <- randomRows(otu_table_gerM, 9)
otu_table_gerF_sub <- randomRows(otu_table_gerF, 9)

#bind together
subset_otu_table <- rbind(otu_table_gerM_sub, otu_table_gerF_sub)
#reset factors
subset_otu_table$var <- factor(subset_otu_table$var)
# run random forest
set.seed(151) 
rf_par <- randomForest( x=subset_otu_table[,1:(ncol(subset_otu_table)-1)] , y=subset_otu_table[ , ncol(subset_otu_table)] , ntree=10000, importance=TRUE, proximities=TRUE )  

# Call:
#  randomForest(x = subset_otu_table[, 1:(ncol(subset_otu_table) -      1)], y = subset_otu_table[, ncol(subset_otu_table)], ntree = 10000,      importance = TRUE, proximities = TRUE)
#                Type of random forest: classification
#                      Number of trees: 10000
# No. of variables tried at each split: 5

#         OOB estimate of  error rate: 33.33%
# Confusion matrix:
#                    sanguisugaNegative sanguisugaPositive class.error
# sanguisugaNegative                  6                  3   0.3333333
# sanguisugaPositive                  3                  6   0.3333333

rf.significance(x=rf_par, xdata=subset_otu_table[,1:(ncol(subset_otu_table)-1)], nperm=1000, ntree=1000)
# Number of permutations:  1000
# p-value:  0.057
# Model not signifiant at p = 0.057
#    Model OOB error:  0.3333333
#    Random OOB error:  0.5555556
#    min random global error: 0.1111111
#    max random global error:  1
#    min random within class error: NA
#    max random within class error:  NA





