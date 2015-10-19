# microarray.r

readInMicroarray <- function(filename){

	microarrayFile <- read.table(filename, header = TRUE)

	M <- microarrayFile[,-dim(microarrayFile)[2]]
	ids <- microarrayFile[,dim(microarrayFile)[2]]
	rownames(M) <- ids
	return(M)
}

getDays <- function(M) {
	splitHeads <- strsplit(colnames(M), split="[.]")
	days <- unique( sapply(splitHeads,head,n=1))
	dayNums <- sapply(sapply( strsplit(days,"[d]"), tail, n=1), as.numeric)
	days <- days[order(dayNums)]
	return(days)
}

getRats <- function(M) {
	splitHeads <- strsplit(colnames(M), split="[.]")
	Rats <- unique( sapply(splitHeads,tail,n=1))
	ratNums <- sapply(sapply( strsplit(Rats,"[R]"), tail, n=1), as.numeric)
	Rats <- Rats[order(ratNums)]
	return(Rats)
}

getRats <- function(M) {
	splitHeads <- strsplit(colnames(M), split="[.]")
	Rats <- unique( sapply(splitHeads,tail,n=1))
	ratNums <- sapply(sapply( strsplit(Rats,"[R]"), tail, n=1), as.numeric)
	Rats <- Rats[order(ratNums)]
	return(Rats)
}



splitMicroarrayMatrixByDay <- function(M) {
	days <- getDays(M)

	colsByDay <- list(mode="matrix")
	for (aDay in days){
		colsByDay[[aDay]] <- M[,sapply(strsplit(colnames(M),"[.]"),head,n=1) %in% c(aDay)]
		cn <- colnames( colsByDay[[aDay]])
		cnNum <- as.numeric( sapply(strsplit( sapply( strsplit( cn, "[.]"), tail,n=1), "[R]"), tail, n=1))
		colsByDay[[aDay]] <- colsByDay[[aDay]][,order(cnNum)]
	}
	return( colsByDay)
}

splitMicroarrayMatrixByRat <- function(M) {
	rats <- getRats(M)

	colsByRat <- list(mode="matrix")
	for (aRat in rats){
		if( aRat != "R5"){
			colsByRat[[aRat]] <- M[,sapply(strsplit(colnames(M),"[.]"),tail,n=1) %in% c(aRat)]
			cn <- colnames( colsByRat[[aRat]])
			cnNum <- as.numeric( sapply(strsplit( sapply( strsplit( cn, "[.]"), head,n=1), "[d]"), tail, n=1))
			colsByRat[[aRat]] <- colsByRat[[aRat]][,order(cnNum)]
		}
	}
	return( colsByRat)
}

loadProbeToGeneMap <- function(filename){
	raw <- read.table(filename, header = TRUE)
	probes <- raw[,1]
	map <- raw[,2,drop=FALSE]
	rownames(map) <- probes
	return(map)
}

probeToGene <- function(probe, pToG){
	return( pToG[probe,])
}

convertProbeDataToGeneData <- function(P, pToG){
	cp <- NCOL(P)
	numcol <- cp + 1
	G <- matrix(0.0,0,numcol)
	colnames(G) <- c(colnames(P), "num_probes")
	mP <- cbind(rownames(P),P)

	rowfunc <- function(r){

		probe <- r[1]
		gene <- as.character( pToG[probe,])
		pexprow <- as.numeric( r[2:(cp+1)])

		if( gene %in% rownames(G)){
			generow <- G[gene,]

			n <- generow[numcol]
			exprow <- generow[-numcol]
			exprow <- (n * exprow) + pexprow
			exprow <- exprow / (n+1)
			generow <- c(exprow, n+1)
			G[gene,] <<- generow

		} else {
			exprow <- pexprow
			generow <- matrix(c(exprow,1),nrow=1)
			rownames(generow) <- gene
			G <<- rbind(G,generow)
		}
	}



	apply(mP,1,rowfunc)

	return(G)
}

