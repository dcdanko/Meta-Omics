
getDayZeroVariableGenes <- function(ngenes=1000){
	getMicroData()
	d <- splitMicroarrayMatrixByDay(g)
	delta <- getTopVaryingGenesComplex(ngenes)

	d0 <<- d$d0[rownames(delta),]

	return(d0)
}

getEveryDayVariableGenes <- function(ngenes=1000){
	getMicroData()
	days <<- splitMicroarrayMatrixByDay(g)

	delta <- getTopVaryingGenesComplex(ngenes)

	days$d0 <<- days$d0[rownames(delta),]
	days$d1 <<- days$d1[rownames(delta),]
	days$d3 <<- days$d3[rownames(delta),]
	days$d6 <<- days$d6[rownames(delta),]
	days$d14 <<- days$d14[rownames(delta),]
	days$d21 <<- days$d21[rownames(delta),]
	days$d28 <<- days$d28[rownames(delta),]
}
getEveryStageVariableGenes <- function(ngenes=1000){
	if(!exists('days')){
		getEveryDayVariableGenes(ngenes=ngenes)
	}

	days$rise <<- cbind(days$d0,days$d1)
	gdays$peak <<- cbind(days$d6,days$d14)
	days$all <<- cbind(days$d0,days$d1,days$d3,days$d6,days$d14,days$d21,days$d28)

}

getRatOneVariableGenes <- function(ngenes=1000){
	getMicroData()
	r1 <- splitMicroarrayMatrixByRat(g)$R1
	delta <- getTopVaryingGenesComplex(ngenes)

	rat1 <<- r1[rownames(delta),]
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

getRatAveVariableGenes <- function(ngenes=1000){
	getMicroData()
	r <- splitMicroarrayMatrixByRat(g)
	rA <- r$R1+r$R2+r$R3+r$R4
	rA <- rA / 4

	delta <- getTopVaryingGenesComplex(ngenes)

	ratA <<- rA[rownames(delta),]
	return(ratA)
}

getRatAveVariableGenesWithInferredDays <- function(ngenes=1000){
	if(!exists('ratA')){
		getRatAveVariableGenes(ngenes)
	}
	source('time-series-clusterer.r')
	inf.ratA <<- inferPoints(ratA)
	return(inf.ratA)
}

getMicroData <- function(){
	if(!exi	ratA <<- rA[rownames(delta),]
sts('g')){
		source('microarray.r')
		print("Reading in gene microarray")
		g <<- readInMicroarray('data/sample_gene_profile.matrix')
	}
}

getProbeMicroData <- function(){
	if(!exists('p')){
		source('microarray.r')
		print("Reading in probe microarray")
		p <<- readInMicroarray('data/sample_probe_profile.matrix')
	}
	if(!exists('pToG')){
		source('microarray.r')
		print('Reading in probe to gene map')
		pToG <<- loadProbeToGeneMap('data/probe2gene.map')
	}
}

getTopVaryingGenes <- function(ngenes=1000){
	getMicroData()
	r <- splitMicroarrayMatrixByRat(g)
	rA <- r$R1+r$R2+r$R3+r$R4
	rA <- rA / 4

	delta <- abs(rA[,1,drop=FALSE] - rA[,5,drop=FALSE] )
	delta <- delta[ order(delta[,1]),,drop=FALSE]
	delta <- delta[ dim(delta)[1]:1,,drop=FALSE]
	delta <- delta[1:ngenes,,drop=FALSE]
	return(delta)

}

getTopVaryingGenesComplex <- function(ngenes=1000){
	getMicroData()
	r <- splitMicroarrayMatrixByRat(g)
	rA <- r$R1+r$R2+r$R3+r$R4
	rA <- rA / 4
	rA.d0 <- rA[,1,drop=FALSE]
	baseline <- cbind(rA.d0,rA.d0,rA.d0,rA.d0,rA.d0,rA.d0,rA.d0)
	rA.logfold <- log(rA / baseline)

	rowf <- function(row){
		return( sum(abs(row)))
	}

	deltavec <- apply(rA.logfold,1,rowf)
	delta <- matrix(deltavec,ncol=1)
	rownames(delta) <- names(deltavec)
	delta <- delta[ order(delta[,1]),,drop=FALSE]
	delta <- delta[ dim(delta)[1]:1,,drop=FALSE]
	delta <- delta[1:ngenes,,drop=FALSE]
	return(delta)

}

getRatAveFoldChange <- function(ngenes=1000){
	if(!exists('ratA')){
		getRatAveVariableGenes(ngenes)
	}

	rA.d0 <- ratA[,1,drop=FALSE]
	baseline <- cbind(rA.d0,rA.d0,rA.d0,rA.d0,rA.d0,rA.d0,rA.d0)
	ratA.fold <<- ratA / baseline
}




