'''Merge two text files with sample names as columns and row names as taxonomy strings'''
import pandas as pd 

dat1 = pd.read_csv("../sequence_taxonomy_table.16s.merged.txt", sep="\t")
dat1 = dat1.groupby("taxonomy").sum()
dat2 = pd.read_csv("../johnston2016_topsoil/sequence_taxonomy_table.16s.filtered.txt", sep="\t")
dat2 = dat2.groupby("taxonomy").sum()
dat3 = pd.read_csv("../ross2018_mammalian_skin/sequence_taxonomy_table.16s.filtered.txt", sep="\t")
dat3 = dat3.groupby("taxonomy").sum()
merged = pd.DataFrame.merge(dat1, dat2, how="outer", left_on="taxonomy", right_on="taxonomy")
merged2 = pd.DataFrame.merge(merged, dat3, how="outer", left_on="taxonomy", right_on="taxonomy")
merged2 = merged2.fillna(0.0)
merged2 = merged2.drop_duplicates()

#write to file
with open("collapsed_otu_table.txt", "w") as outfile: 
		merged2.to_csv(outfile, sep="\t")