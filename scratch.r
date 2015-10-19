
# source("microarray.r")
# pg <- loadProbeToGeneMap("data/probe2gene.map")
# m <- readInMicroarray("data/sample_probe_profile.matrix")
# r <- splitMicroarrayMatrixByRat(m)
# r1 <- r$R1
# r1.gene <- convertProbeDataToGeneData(r1,pg)
# r1.gene <- r1.gene[,1:7]
# r1.base <- r1.gene[,1]
# r1.base.m <- cbind(r1.base,r1.base,r1.base)
# r1.base.m <- cbind(r1.base.m,r1.base.m,r1.base)
# r1.delta <- r1.gene - r1.base

# data <- t(r1.delta)
# data <- t(r1.delta[1:5000,])
# source("correlation-matrix.r")

rowinterp <- function(r){
  #                   d0   d1   d2 d3   d4 d5 d6  
  infr <- na.spline(c(r[1],r[2],NA,r[3],NA,NA,r[4],
  #                   d7 d8 d9 d10 d11 d12 d13 d14
                      NA,NA,NA,NA, NA, NA, NA, r[5],
  #                   d15 d16 d17 d18 d19 d20 d22
                      NA, NA, NA, NA, NA, NA, r[6],
  #                   d23 d24 d25 d26 d27 d28                    
                      NA, NA, NA, NA, NA, r[7]))
  return(infr)
}

 c <- matrix(apply(r1.delta[1:2,],1,FUN=rowinterp),byrow=TRUE,nrow=2)


> source("clustering.r")
> d.r1.infdel <- makeDistMatrix(r1.delta.inf
r1.delta.inf
> d.r1.infdel <- makeDistMatrix(r1.delta.inf, 'dtw')
^C
^C  
^C^C> 
> d.r1.infdel <- makeDistMatrix(r1.delta.inf[1:100,], 'dtw')
> d.r1.5 <- kmeans(d.r1.infdel,5)
> clustering <- d.r1.5
> series <-r1.
-r1.base       -r1.base.m     -r1.delta      -r1.delta.inf  -r1.gene
> series <- r1.delta.inf
> source('workflows.r')
Error: object 'cnum' not found
> cnum <- -1
> source('workflows.r')
Error in colors[iclust] : invalid subscript type 'list'
> colors
[1] "#FF0000FF" "#CCFF00FF" "#00FF66FF" "#0066FFFF" "#CC00FFFF"
> colors[2]
[1] "#CCFF00FF"
> nseries
[1] 17330
> clustering <- d.r1.5[1:100,]
Error in d.r1.5[1:100, ] : incorrect number of dimensions
> clustering <- d.r1.5
> series <- r1.delta.inf[1:100,]
> d.r1.5$cluster[3]
0610005H09RIK 
            2 
> clustering <- d.r1.5$cluster
> source('workflows.r')
> 