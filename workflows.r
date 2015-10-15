# workflows.r

# from data file to distance matrix

source("~/Meta-Omics/microarray.r")
source("~/Meta-Omics/clustering.r")
library("dtw")

m <- readInMicroarray("~/Dropbox/Microarray_data/sample_probe_profile.matrix")
r <- splitMicroarrayMatrixByRat(m)

r <- r$R1
r.c <- r - apply(r,1, function(row) median(range(row))[1])
numrowstouse <- 100
rise <- r[1:numrowstouse,1:5]
# d.c <- makeDistMatrix(r.c[1:1000,], "dtw")

# # Do K-Means Clustering

# d.c.10 <- kmeans(d.c,10)


# plot a set of time-series as lines
# start <- 1
# end <- 1000
# clustering <- rise.h.100
# series <- rise[start:end,]

# time <- c(0,1,3,6,14) #,21,28)

# nclust <- 100 # dim(clustering$centers)[1]
# nseries <- dim(series)[1]

# plot(range(time), range(series), type='n', xlab="Day", ylab="Exp. Level")
# colors <- rainbow( nclust)

# for (i in 1:nseries){
# 	iclust <- clustering[i] # clustering$cluster[i]

# 	if(cnum == -1 || cnum == iclust ){
# 		lines(time, series[i,], col=colors[iclust])
# 	}
# }

# title("Expression level by day")

# Center rows

# m.centered <- m - apply(m,1, function(r) median(range(r))[1])

# # scale rows

# m.scaled <- m / apply(m,1,max)

# library("pvclust")

# p.rise <- pvclust(t(rise), method.hclust="ward.D2", method.dist="dtw")

# plot(p.rise)
# pvrect(p.rise, alpha=0.80, color="green")

