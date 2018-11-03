#zach = read.csv("http://www-personal.umich.edu/~mejn/dcsbm/ZacharyCorrectOutput/DegreeCorrected/EdgeLists.tsv", sep = "\t")
#if multiedges
#multiedge$count = 1
#aggregate(multiedge$count, by = list(multiedge$X0, multiedge$X1), FUN = "sum")

library(ergm.count)
data(zach)
zach = data.frame(as.edgelist(zach))
save(zach, file="sbm/data/zach.rda")
zach$count = (1:13) #recycle
zach_weighted = zach
save(zach_weighted, file="sbm/data/zach_weighted.rda")
