.summary.SeqCNAInfo = function(object) {
	cat("Basic information:\n")
	cat(paste("  SeqCNAInfo object with ", length(object@tumour[[1]]), " ",object@win,"Kbp-long windows.\n", sep=""))
	has.pem = (length(object@normal)==0 && length(object@tumour)>1) || length(object@normal)>1
	cat(paste("  PEM information is ", ifelse(has.pem,"","not "), "available.\n", sep=""))
	cat(paste("  Paired normal is ", ifelse(length(object@normal)>0,"","not "), "available.\n", sep=""))
	seqtab = rle(object@seq)$values
	if (is.null(object@build)) {
		build = "unknown"
	} else {
		if (is.na(object@build) || nchar(object@build)==0)
			build = "unknown"
		else
			build = object@build
	}
	cat(paste("  Genome and build ",build,
		" (chromosomes ", head(seqtab,1), " to ", tail(seqtab,1),").\n", sep=""))
	if (length(object@skip) > 0)
		cat(paste("Total filtered windows: ", length(object@skip), ".\n", sep=""))
	cat(paste("The profile is ", ifelse(ncol(object@output)>0,"","not yet "), "normalized and ",
		ifelse(ncol(object@output)>3,"","not yet "), "segmented.\n", sep=""))
	if (length(object@minCN) > 0) {
		mincn = object@minCN
		maxcn = mincn + length(object@thr)
		cat(paste("Copy numbers called from ", mincn, " to ", maxcn,".\n", sep=""))
		for (cn in mincn:maxcn)
			cat(paste("  CN", cn, ":", length(which(object@output$CN==cn))," windows.\n", sep=""))
	}
}
