#should be done with a binary grouping

import pandas as pd
from skbio.stats.composition import ancom


dat = pd.DataFrame.from_csv("sequence_table.16s.filtered.txt", sep="\t", index_col=0, header=0)
clust = pd.DataFrame.from_csv("philr_cluster.txt", sep="\t", index_col=0, header=0)
meta = pd.DataFrame.from_csv("test.txt", sep="\t", index_col=0, header=0)


# merge = pd.DataFrame.merge(dat, clust, left_index=True, right_index=True)
# a = list(merge['cluster'])
# b = list(merge['states'])
# grouping = pd.Series(a, index=b)

merge = pd.DataFrame.merge(dat, meta, left_index=True, right_index=True)
a = list(merge['MachRes1'])
grouping = pd.Series(a, index=merge.index.values)

#pseudocount
numeric_cols = [col for col in dat if dat[col].dtype.kind !='0']
dat[numeric_cols] += 1

#ancom analysis
results = ancom(dat, grouping)
print(results)