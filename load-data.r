
getDayZeroVariableGenes <- function(ngenes=1000){
	getMicroData()
	d <- splitMicroarrayMatrixByDay(m)
	delta <- getTopVaryingGenes(ngenes)

	d0 <- convertProbeDataToGeneData(d$d0[rownames(delta),], pToG)
	d0 <<- d0[,-ncol(d0)]

	# r1 <- convertProbeDataToGeneData(r$R1[rownames(delta),1,drop=FALSE] ,pToG)
	# print("Finished converting 1/5 probe data.")
	# r2 <- convertProbeDataToGeneData(r$R2[rownames(delta),1,drop=FALSE] ,pToG)
	# print("Finished converting 2/5 probe data.")
	# r3 <- convertProbeDataToGeneData(r$R3[rownames(delta),1,drop=FALSE] ,pToG)
	# print("Finished converting 3/5 probe data.")
	# r4 <- convertProbeDataToGeneData(r$R4[rownames(delta),1,drop=FALSE] ,pToG)
	# print("Finished converting 4/5 probe data.")
	# r5 <- convertProbeDataToGeneData(m[rownames(delta),7,drop=FALSE] ,pToG)
	# print("Finished converting 5/5 probe data.")

	# r1 <- r1[,-NCOL(r1)]
	# r2 <- r2[,-NCOL(r2)]
	# r3 <- r3[,-NCOL(r3)]
	# r4 <- r4[,-NCOL(r4)]
	# r5 <- r5[,-NCOL(r5)]

	# d0 <<- cbind(r1, r2, r3, r4, r5)
	return(d0)
}

getRatOneVariableGenes <- function(ngenes=1000){
	getMicroData()
	r1 <- splitMicroarrayMatrixByRat(m)$R1
	delta <- getTopVaryingGenes(ngenes)

	rat1 <- convertProbeDataToGeneData(r1[rownames(delta),] ,pToG)
	rat1 <<- rat1[,-ncol(rat1)]
	return(rat1)
}

getRatOneVariableGenesWithInferredDays <- function(ngenes=1000){
	if(!exists('rat1')){
		getRatOneVariableGenes(ngenes)
	}
	source('time-series-clusterer.r')
	inf.rat1 <<- inferPoints(rat1)
	return(inf.rat1)
}

getMicroData <- function(){
	if(!exists('m')){
		source('microarray.r')
		print("Reading in microarray")
		m <<- readInMicroarray('data/sample_probe_profile.matrix')
	}
	if(!exists('pToG')){
		source('microarray.r')
		print('Reading in probe to gene map')
		pToG <<- loadProbeToGeneMap('data/probe2gene.map')
	}
}

getTopVaryingGenes <- function(ngenes=1000){
	getMicroData()
	r <- splitMicroarrayMatrixByRat(m)
	r1 <- r$R1
	rA <- r$R1+r$R2+r$R3+r$R4
	rA <- rA / 4

	delta <- abs(rA[,1,drop=FALSE] - rA[,5,drop=FALSE] )
	delta <- delta[ order(delta[,1]),,drop=FALSE]
	delta <- delta[ dim(delta)[1]:1,,drop=FALSE]
	delta <- delta[1:ngenes,,drop=FALSE]
	return(delta)

}




