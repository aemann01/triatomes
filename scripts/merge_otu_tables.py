'''Merge two text files with sample names as columns and row names as taxonomy strings'''
import pandas as pd 

dat1 = pd.read_csv("sourcetracker/sequence_table_for_sourcetracker.txt", sep="\t")
dat2 = pd.read_csv("sequence_taxonomy_table.16s.merged.txt", sep="\t")
merged = pd.DataFrame.merge(dat1, dat2, how="outer", left_on="Taxonomy", right_on="taxonomy", indicator=True)

#now collapse down based on row names and get sum
collapse_data = merged.groupby("Taxonomy").sum()

#write to file
with open("sourcetracker/collapsed_otu_table.txt", "w") as outfile: 
		collapse_data.to_csv(outfile, sep="\t")