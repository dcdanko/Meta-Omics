
getDayZeroVariableGenes <- function(ngenes=1000){
	getMicroData()
	d <- splitMicroarrayMatrixByDay(m)
	delta <- getTopVaryingGenes(ngenes)

	d0 <- convertProbeDataToGeneData(d$d0[rownames(delta),], pToG)
	d0 <<- d0[,-ncol(d0)]

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




