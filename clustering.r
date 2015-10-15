# clustering.r

makeDistMatrixCustom <- function(data, dFunc) {
	nrows <- dim(data)[1]

	d <- matrix(0,nrow=nrows,ncol=nrows)
	for (i in 1:nrows){
		iVec <- data[i,]
		for (j in 1:i){
			if( j!=i){
				jVec <- data[j,]
				dist <- dFunc(iVec,jVec)
			} else {
				dist <- 0
			}
			d[i,j] <- dist
			d[j,i] <- dist
		}
	}
	colnames(d) <- rownames(data)
	rownames(d) <- rownames(data)
	return(as.dist(d))
}

makeDistMatrix <- function(data, meth) {
	d <- as.dist( dist(data, method=meth))
	return(d)
}


manhattanDist <- function(v1,v2) {
	return( sum(abs(v1-v2)))
}

euclideanDist <- function(v1,v2) {
	return( sqrt( sum( (v1-v2)*(v1-v2))))
}


centerAndScaleRows <- function(m){
	m.centered <- m - apply(m,1, function(r) median(range(r))[1])
	m.scaled <- m.centered / apply(m.centered,1,max)
	return(m.scaled)
}
