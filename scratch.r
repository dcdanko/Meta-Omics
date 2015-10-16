> setwd("Meta-Omics/")
> source("microarray.r")
> m <- readInMicroarray("data/sample_probe_profile.matrix")
> pg <- loadProbeToGeneMap("data/probe2gene.map")
> 
> r <- splitMicroarrayMatrixByRat(m)
> rA <- r$R1+r$R2+r$R3+r$R4
> rA <- rA / 4
> delta <- abs(rA[1]-rA[5])
> delta <- delta[ order(delta[,1]),,drop=FALSE]
> delta <- delta[ dim(delta)[1]:1,,drop=FALSE]
k <- delta[1:1000,,drop=FALSE]
> k5 <- m[,7,drop=FALSE][rownames(k),,drop=FALSE]
> k1 <- r$R1[rownames(k),1,drop=FALSE]
> k2 <- r$R2[rownames(k),1,drop=FALSE]
> k3 <- r$R3[rownames(k),1,drop=FALSE]
> k4 <- r$R4[rownames(k),1,drop=FALSE]
> k.d0 <- cbind(k1,k2,k3,k4,k5)
> k.d0.gene <- convertProbeDataToGeneData(k.d0,pg)
> k.d0.gene <- k.d0.gene[,1:5]
> data <- t(k.d0.gene)
> source("correlation-matrix.r")
> anet <- graph.adjacency(abs(adja),mode="undirected",weighted=TRUE,diag=FALSE)
> catch <- page.rank(anet,vids=V(anet))
> top <- catch$vector[order(catch$vector)][850:1]
> top[1:50]
           MMP3            DIO1           H2-M2          CRELD2       LOC195359 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
           OSTB            CTSZ            MMP7         PPP1R9A          H2-T23 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
        TNFAIP2   6330415F13RIK   9130208D14RIK         UGT1A10             VIP 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
         PARP14           ACTN4            LCP2       LOC383540             MIF 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
         PHLDA1           STIP1   2310075E07RIK            SDC2             ELN 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
           RERG          MS4A6B           H2-Q2   1300018P11RIK SCL0001259.1_60 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
  E030026I10RIK   2310046K10RIK             EMB   0610010I05RIK        AI481105 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
          CHODL          MAPK13            THBD           PLEK2          DEFB27 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001746812 
          PRKCZ           MYO1A   5730577G12RIK           CIDEB           CXCL5 
    0.001746812     0.001746812     0.001746812     0.001746812     0.001625809 
       SLC30A10         CYP2C55         FPR-RS2   2610041P08RIK          SAMHD1 
    0.001625809     0.001625809     0.001625809     0.001625809     0.001625809 
> ns <- names(top[1:20])
> vs <- V(anet)[name %in% ns]
> sub <- induced.subgraph(anet,vs)
> plot(sub)

