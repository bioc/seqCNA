##################################################################
### MAIN FUNCTIONS
##################################################################
runSeqsumm = function(summ.win=50, file=NULL, folder=NULL, output.file="seqsumm_out.txt", samtools.path="samtools") {
	if (! is.null(file)) {
		if (! file.exists(file)) {
			stop("File does not exist.")
		} else {
			folder = dirname(file)
			filename = basename(file)
			file.type = toupper(gsub(".*[.]","",filename))
			if (! file.type%in%c("SAM","BAM"))
				stop("File should have either SAM or BAM extension.")
		}
	} else {
		if (! .folderExists(folder)) 
			stop("Folder does not exist.")
	}
	wd = getwd()
	setwd(folder)
	if (length(list.files(pattern="^seqsumm_out.txt$")) > 0)
		warning("A seqsumm output file already exists in the folder and will be overwritten.", immediate.=TRUE)
	if (is.null(file)) {
		filename = list.files(pattern=".sam")[1]
		if (is.na(filename)) {
			file.type = "BAM"
			filename = list.files(pattern=".bam")[1]
			if (is.na(filename)) {
				setwd(wd)
				stop("There are not any SAM or BAM files in the specified folder.")
			}
		} else {
			file.type = "SAM"
		}
	}
	cmd = ifelse(file.type=="BAM", paste(samtools.path," view ",filename,sep=""), paste("cat",filename))
	message(paste("Reading records from", filename))
	.C('seqsumm', as.integer(summ.win*1000), cmd, output.file)
	setwd(wd)
}

readSeqsumm = function(build="", tumour.data=NULL, normal.data=NULL, folder=NULL, normal.folder=NULL, resample.win=NULL, sex=TRUE, nproc=2) {
	.readSeqsumm = function(tab, columns=c(5:ncol(tab))) {
		if (ncol(tab) < 5)
			stop("The seqsumm file has a wrong format or is corrupt.")
		rco = list()
		m = min(columns)
		for (chr in chrs) {
			chrtab = tab[tab$chrom==chr,]
			for (i in columns) {
				chrtabti = chrtab[,i]
				chrtabti = chrtabti[1:lengths[chr]]  # make them of same length as the chromosome
				chrtabti[is.na(chrtabti)] = 0  # fill sample windows with 0s if necessary
				if (length(rco) < i-m+1)
					rco[[i-m+1]] = chrtabti
				else
					rco[[i-m+1]] = c(rco[[i-m+1]], chrtabti)
			}
		}
		rco
	}
	# check if build is present and supported by annotation
	if (length(build)==0)
		build = ""
	supported = FALSE
	if (! is.null(build))
		if (! is.na(build) && nchar(build)>0)
			if (build %in% supported.builds())
				supported = TRUE
	if (! supported)
		message("Note: chromosome lengths will be estimated from sample.")
	message("Reading summarized data...")
	flush.console()
	# read tumoural data
	if (! is.null(tumour.data)) {
		tab.tumour = as.data.frame(tumour.data)
	} else {
		if (.folderExists(folder)) {
			if (.outputExists(folder))
				tab.tumour = read.table(file.path(folder, "seqsumm_out.txt"), header=TRUE, comment.char="")
			else
				stop("No seqsumm output exists in the folder for the tumoural sample.")
		} else {
			stop("The folder for the tumoural sample does not exist.")
		}
	}
	# read normal data if applicable
	if (! is.null(normal.data)) {
		tab.normal = as.data.frame(normal.data)
	} else {
		tab.normal = NULL
		if (length(normal.folder) > 0 && !.folderExists(normal.folder))
			warning("The folder for the normal sample does not exist. Only the tumoural sample will be read.")
		if (.folderExists(normal.folder)) 
			if (.outputExists(normal.folder))
				tab.normal = read.table(file.path(normal.folder, "seqsumm_out.txt"), header=TRUE, comment.char="")
	}
	# perform downsampling if requested
	if (!is.null(resample.win)) {
		message("Resampling summarized data...")
		flush.console()
		tab.tumour = .resample(tab.tumour, resample.win, nproc)
		if(!is.null(tab.normal))
			tab.normal = .resample(tab.normal, resample.win, nproc)
	}
	# calculate window size
	win = diff(tab.tumour$win.start[1:2])/1000
	if (win == 0)
		stop("Window size is zero.")
	tab.tumour$chrom = as.character(tab.tumour$chrom)
	chrs = as.character(rle(tab.tumour$chrom)$values)
	# remove sex chromosomes if requested
	if (! sex) {
		if (build %in% c("hg17","hg18","hg19","hg20"))
			sex.chrs = c("X","chrX","Y","chrY","23","chr23","24","chr24")
		else
			sex.chrs = c("X","chrX","Y","chrY")
		wh = which(chrs%in%sex.chrs)
		if (length(wh) > 0) {
			chrs = chrs[-wh]
			tab.tumour = tab.tumour[-which(tab.tumour$chrom%in%sex.chrs),]
			if (! is.null(tab.normal))
				tab.normal = tab.normal[-which(tab.normal$chrom%in%sex.chrs),]
		}
	}
	# generate chromosome window indices
	if (supported) {
		info = paste(build,"_len",sep="")
		data(list=info)
		info.chrs = as.character(as.vector(get(info)[,1]))
		info.lengths = as.integer(as.vector(get(info)[,2]))
	} else {
		info.chrs = chrs
		info.lengths = sapply(info.chrs, function(i) max(tab.tumour$win.start[tab.tumour$chrom==i])+win*1000)
	}
	lengths = sapply(chrs, function(i) ceiling(info.lengths[info.chrs==i]/(win*1000)))
	names(lengths) = chrs
	all.chrs = rep(chrs, lengths)
	all.pos = unlist(lapply(chrs, function(i) c(1:lengths[i])))
	# create SeqCNAInfo object from input and calculated information
	rcount.tumour = .readSeqsumm(tab.tumour)
	if (!is.null(tab.normal)) {
		rcount.normal = .readSeqsumm(tab.normal)
		reads.gc = .readSeqsumm(tab.normal,3)[[1]]
		reads.mapq = .readSeqsumm(tab.normal,4)[[1]]
	} else {
		rcount.normal = list()
		reads.gc = .readSeqsumm(tab.tumour,3)[[1]]
		reads.mapq = .readSeqsumm(tab.tumour,4)[[1]]
	}
	rco = new("SeqCNAInfo", tumour=rcount.tumour, normal=rcount.normal, gc=reads.gc, mapq=reads.mapq,
		seq=all.chrs, pos=all.pos, build=build, win=win)
}

applyFilters = function(rco, pem.filter=0, trim.filter=0, mapp.filter=0, mapq.filter=0, cnv.filter=FALSE, plots=TRUE, folder=NULL, nproc=2) {
	stop.text = "Run readSeqsumm first!"
	if (length(rco@tumour) == 0) stop(stop.text)
	if (length(rco@seq) == 0) stop(stop.text)
	if (length(rco@pos) == 0) stop(stop.text)
	if (length(rco@win) == 0) stop(stop.text)
	use.normal = length(rco@normal)!=0
	chrs = rle(rco@seq)$values
	if (is.null(rco@build)) {
		supported = FALSE
	} else {
		if (is.na(rco@build) || nchar(rco@build)==0)
			supported = FALSE
		else
			supported = rco@build%in%supported.builds()
	}
	which.cnv = which.mapp = which.mapq = c()
	which.filter = which.imp = which.q = which.trim = c()
	if (supported) {
		if (!use.normal || mapp.filter>0 || cnv.filter) {
			message("Loading genome build information...")
			data(list=rco@build)
			data(list=paste(rco@build,"_len",sep=""))
			genome.data = .buildGenomeInfo(get(rco@build), get(paste(rco@build,"_len",sep="")), chrs, rco@win, nproc=nproc)
		}
		if (use.normal) {
			x = rco@normal[[1]]
		} else {
			x = c()
			for (chr in chrs)
				x = c(x, as.numeric(as.vector(genome.data$GC[which(genome.data$Chr==chr)])))
			x[x==-1] = NA # if necessary, GC windows with code -1 recoded to NAs
		}
		# FILTER WINDOWS WITH CNVs or low mappability
		if (cnv.filter) {
			if (use.normal) {
				cnv.filter = FALSE
				message("Note: CNV filter not meaningful with paired-normal. Disabled.")
			} else {
				which.cnv = which(genome.data$CNV>0.95)
			}
		}
		if (mapp.filter)
			which.mapp = which(genome.data$Mapp<mapp.filter)
	} else { # avoid CNV filter, and estimate GC content and mappability
		if (use.normal) {
			x = rco@normal[[1]]
		} else {
			x = rco@gc
			x[x<=0] = NA
			x[x>=1] = NA
			message("Note: GC content will be estimated from sample.")
		}
		if (cnv.filter) {
			message("Note: the genome build does not support CNV-based filtering.")
			cnv.filter = FALSE
		}
		if (mapp.filter) {
			message("Note: the genome build does not support mappability-based filtering.")
			mapp.filter = 0
		}
	}
	y = rco@tumour[[1]]

	# (PRE-)NORMALIZE THE PROFILE FOR BETTER FILTER ASSESSMENT AND IMPROVED TRIM-BASED FILTERING
	if (! use.normal) { # pre-normalize tumoural against GC
		pre.x = x
		pre.y = y
	} else { # or normal against GC
		pre.x = rco@gc
		pre.y = x
	}
	method = "quadratic"
	reg.method = ifelse(method=="loess", "loess", "lm")
	which.ok = intersect(intersect(which(pre.y>0),which(!is.na(pre.y))),which(!is.na(pre.x)))
	x.pos = pre.x[which.ok]
	y.pos = pre.y[which.ok]
	lomo = eval(call(reg.method, .makeFormula("x.pos","y.pos",method,use.normal,environment())))
	y.norm = rep(NA, length(y))
	z.pos = lomo$residuals / lomo$fitted
	y.norm[which.ok] = z.pos
	# FILTER WINDOWS WITH LOW MEAN MAPPING QUALITY
	if (mapq.filter)
		which.mapq = intersect(which(rco@mapq<mapq.filter), which.ok)
	which.filter = unique(c(which.mapq, which.cnv, which.mapp))
	# FILTER WINDOWS WITH HIGH RATIO OF IMPROPER PAIR-END MAPPINGS IN THE NORMAL (TUMOURAL) SAMPLE,
	# A SIGN OF REPETITIVE SEQUENCES, MOSTLY HAPPENING AROUND THE CENTROMERIC REGIONS
	pem.allowed = (use.normal && length(rco@normal) > 1) || length(rco@tumour) > 1
	if (pem.filter && !pem.allowed)
		message("Note: PEM filter will not be applied due to lack of PEM information.")
	message("Applying filters...")
	do.pem.filter = pem.filter && pem.allowed
	if (do.pem.filter) {
		types = c(2,4,5) # type 3 reads are non-existent in some alignments, so better not use them 
		if (use.normal)
			improper = Reduce('+',(lapply(types,function(i)rco@normal[[i]])))/rco@normal[[1]]
		else
			improper = Reduce('+',(lapply(types,function(i)rco@tumour[[i]])))/rco@tumour[[1]]
		thr = quantile(improper[is.finite(improper)], probs=1-pem.filter)
		which.imp = which(improper > thr)
		which.filter = union(which.filter, which.imp)
	} else {
		pem.filter = 0
	}
	# FILTER WINDOWS WITH EXTREME COUNTS, PRONE TO RESULT IN GREATER NOISE IN THE RATIO PROFILE
	if (trim.filter[1]) {
		if (use.normal && trim.filter[1]==1) {
			trim.filter[1] = 0
			message("Automatic trimming is only available without normal sample. Disabled.")
		}
		trim.val = ifelse (!use.normal&&trim.filter[1]==1, .autoTrim(z.pos,nproc=nproc), 1-trim.filter[1])
		trim.val[2] = ifelse (length(trim.filter) > 1, trim.filter[2], 1-trim.val[1])
		quants = quantile(z.pos, rev(trim.val))
		which.trim = which(y.norm<quants[1] | y.norm>quants[2])
		if (length(which.ok)>0)
			which.trim = union(which.trim, c(1:length(y))[-which.ok])
		which.q = intersect(which.trim, which.ok)
		which.filter = union(which.filter, which.trim)
	}
	which.nok = c(1:length(y))[-which.ok]
	# OUTPUT FILTERING RESULTS
	zeroes.len = length(y)-length(z.pos)
	message("  Windows without reads or info: ",appendLF=FALSE)
	message(length(which.nok))
	message("  Total filtered windows:        ",appendLF=FALSE)
	message(length(union(which.filter,which.nok)))
	message("  Remaining windows:             ",appendLF=FALSE)
	message(length(y)-length(union(which.filter,which.nok)))
	flush.console()
	if (plots) {
		mask  = c(trim.filter[1]>0, mapp.filter>0, mapq.filter>0, pem.filter>0, cnv.filter)
		l = list(which.trim, which.mapp, which.mapq, which.imp, which.cnv)[mask]
		ll = length(unlist(l))
		lf = length(which(mask==TRUE))
		nrows = ifelse(lf>3,2,1)
		ncols = ceiling(lf/nrows)
		if (ll > 0) {
			message("Plotting filtering...")
			flush.console()
			all.cols = c("#9E2862FF","#FF7400FF","#FF00BBFF","#006388FF","#92EC00FF")
			cols = all.cols[mask]
			pchs = c(1:5)[mask]
			if (.folderExists(folder)) {
				jpeg(file.path(folder,"seqCNA_filters.jpg"),1920,1080)
			} else {
				if (lf<4) m = t(matrix(c(1,1,rep(1,ncols),2,2,c(3:(lf+2))),ncol=2))
				else m = t(matrix(c(1,1,1,1,2,2,3,4,2,2,5,6),ncol=3))
				layout(m)
			}
			# profile with filtered windows
			par(oma=c(0,0,0,0))
			old.par = par(mar=c(5,4,4,0))
			if (trim.filter[1]) {
				z = y.norm
				yl = quantile(z.pos[is.finite(z.pos)], c(trim.val[2]/2, 1-(1-trim.val[1])/2))
				if (!use.normal) ylab = "Pre-normalized tumoural profile"
				else ylab = "Normalized normal profile"
			} else {
				z = y
				yl = c(0,quantile(y, 0.9995))
				ylab = "Raw counts"
			}
			binsX = 500
			bandX = ceiling(length(z) / binsX)
			binsY = 200
			bandY = diff(quantile(z[is.finite(z)],c(0.01,0.99))) / binsY
			suppressMessages(suppressWarnings(
			smoothScatter(postPlotHook=NULL, colramp=colorRampPalette(c("white","black")),
				bandwidth=c(bandX,bandY), nbin=c(binsX,binsY),
				z, cex=0.5, pch=19, col="#00000033", ylim=yl, cex.axis=1.5,cex.lab=1.5, ylab=ylab, frame.plot=FALSE)
			))
			for (i in 1:length(l))
				points(l[[i]], z[l[[i]]], cex=1, pch=pchs[i], col=cols[i])
			if (.folderExists(folder))
				dev.off()
			par(mar=old.par)
			# Venn diagrams
			if (lf > 1) {
				venn.mask = c(pem.filter>0, mapp.filter>0, mapq.filter>0, cnv.filter, trim.filter[1]>0)
				venn.l = list(which.imp,intersect(which.mapp,which.ok),intersect(which.mapq,which.ok),which.cnv,which.q)[venn.mask]
				names(venn.l) = c("PEM","mapp.","quality","CNV","trimmed")[venn.mask]
				venn.cols = c("#006388FF","#FF7400FF","#FF00BBFF","#92EC00FF","#9E2862FF")[venn.mask]
				if (.folderExists(folder))
					jpeg(file.path(folder,"seqCNA_filters_overlap.jpg"),1920,1080)
				counts = list(sapply(.overLapper(setlist=venn.l, sep="_", type="vennsets")$Venn_List, length))
				old.par = par(mar=rep(0,4))
				dummy = .vennPlot(counts=counts, lines=venn.cols, lcol=venn.cols, mymain="", mysub="", ccex=1.25)
				par(mar=old.par)
				if (.folderExists(folder))
					dev.off()
				if (.folderExists(folder)) {
					jpeg(file.path(folder,"seqCNA_filters_thresholds.jpg"),1920,1080)
					par(mfrow=c(nrows,ncols))
				}
			}
			# per-filter plots and thresholds
			if (mask[1]) {
				d = density(y.norm[y.norm<quants[2]*1.5 & y.norm>quants[1]*1.5],na.rm=TRUE)
				plot(Inf, xlim=c(min(d$x),max(d$x)), ylim=c(0,max(d$y)),
					xlab=ylab,ylab="Density",main="Trimming filter",cex.lab=1.5,cex.axis=1.5)
				rect(quants[1], -max(d$y)*2, quants[2], max(d$y)*2, col=gsub("FF$","44",all.cols[1]), lty=0)
				abline(h=0, col="grey")
				lines(d, lwd=2)
				abline(v=quants, col=all.cols[1], lwd=2)
				text(quants[1],max(d$y),round(quants[1],2),adj=c(-0.25,0.75), cex=1.5, col=all.cols[1])
				text(quants[2],max(d$y),round(quants[2],2),adj=c(-0.25,0.75), cex=1.5, col=all.cols[1])
			}
			if (mask[2]) {
				d = density(genome.data$Mapp,na.rm=TRUE)
				plot(Inf, xlim=c(0,1), ylim=c(0,max(d$y)),
					xlab="Mappability",ylab="Density",main="Mappability filter",cex.lab=1.5,cex.axis=1.5)
				rect(mapp.filter, -max(d$y)*2, 2, max(d$y)*2, col=gsub("FF$","44",all.cols[2]), lty=0)
				abline(h=0, col="grey")
				lines(d, lwd=2)
				abline(v=mapp.filter, col=all.cols[2], lwd=2)
				text(mapp.filter,max(d$y),mapp.filter,adj=c(-0.25,0.75), cex=1.5, col=all.cols[2])
			}
			if (mask[3]) {
				d = density(rco@mapq,na.rm=TRUE)
				plot(Inf, xlim=c(0,max(d$x)), ylim=c(0,max(d$y)),
					xlab="Mean mapping quality",ylab="Density",main="Mapping quality filter",cex.lab=1.5,cex.axis=1.5)
				rect(mapq.filter, -max(d$y)*2, max(d$x)*2, max(d$y)*2, col=gsub("FF$","44",all.cols[3]), lty=0)
				abline(h=0, col="grey")
				lines(d, lwd=2)
				abline(v=mapq.filter, col=all.cols[3], lwd=2)
				text(mapq.filter,max(d$y),mapq.filter,adj=c(-0.25,0.75), cex=1.5, col=all.cols[3])
			}
			if (mask[4]) {
				d = density(improper[improper<thr*1.5],na.rm=TRUE)
				plot(Inf, xlim=c(0,max(d$x)), ylim=c(0,max(d$y)),
					xlab="Improper/proper read ratio",ylab="Density",main="PEM filter",cex.lab=1.5,cex.axis=1.5)
				rect(-1, -max(d$y)*2, thr, max(d$y)*2, col=gsub("FF$","44",all.cols[4]), lty=0)
				abline(h=0, col="grey")
				lines(d, lwd=2)
				abline(v=thr, col=all.cols[4], lwd=2)
				text(thr,max(d$y),round(thr,2),adj=c(-0.25,0.75), cex=1.5, col=all.cols[4])
			}
			if (.folderExists(folder)) dev.off()
		} else {
			message("There are no filtered windows to plot.")
		}
	}
	rco@y = y
	rco@x = x
	rco@skip = as.numeric(which.filter)
	rco
}

runSeqnorm = function(rco, norm.win=NULL, method="quadratic", lambdabreak=8, minSeg=7, maxSeg=35, nproc=2, plots=TRUE, folder=NULL) {
	if (length(rco@tumour) == 0) stop("Run readSeqsumm first!")
	if (length(rco@y) == 0) stop("Run applyFilters first!")
	if (length(rco@x) == 0) stop("Run applyFilters first!")
	if (!is.null(norm.win)) {
		if ((norm.win/rco@win) %% 1 != 0)
			message("Note: normalization window will be approximated to the closest multiple of the summarization window.")
	} else {
		if (length(rco@y) > 200000) {
			f = round(length(rco@y) / 50000)
			norm.win = rco@win * f
		} else {
			norm.win = rco@win
		}
	}
	rcount.prof = seqnormX(rco@x, rco@y, chr=rco@seq, skip=rco@skip, use.normal=length(rco@normal)!=0,
			resolution=rco@win/norm.win, nproc=nproc,method=method, out.dir=folder, 
			out.jpg=file.path(folder,"seqCNA_normalization.jpg"), plots=plots,
			lambdabreak=lambdabreak, minRegions=minSeg, maxRegions=maxSeg)
	profiles = cbind(rco@seq, (rco@pos-1)*rco@win, round(rcount.prof,4))
	colnames(profiles) = c("chrom", "win.start", "normalized")
	rco@output = data.frame(profiles)
	rco@output$normalized = as.numeric(as.vector(rco@output$normalized))
	rco
}

runGLAD = function(rco, lambdabreak=8, nproc=2) {
	if (length(rco@tumour) == 0) stop("Run readSeqsumm first!")
	if (length(rco@y) == 0) stop("Run applyFilters first!")
	if (ncol(rco@output) == 0) stop("Run runSeqnorm first!")
	message("Running segmentation...")
	flush.console()
	which.ok = which(!is.na(rco@output$normalized))
	chrs = rle(as.vector(rco@output$chrom))$values
	ch = 0
	all.prof = rco@output[, c(3,1)]
	cluster = makeCluster(rep("localhost",nproc), type="SOCK") 
	registerDoSNOW(cluster)
	res = foreach(ch=chrs, .packages="GLAD") %dopar% {
		prof = all.prof[which(all.prof$chrom==ch), ]
		input = cbind(c(1:nrow(prof)),prof)
		colnames(input) = c("PosOrder","LogRatio","Chromosome")
		profileCGH = as.profileCGH(data.frame(input))
		res = glad(profileCGH, verbose=FALSE, bandwidth=1, lambdabreak=lambdabreak)$profileValues$Smoothing
	}
	stopCluster(cluster)
	res = unlist(res)
	rcount.prof.segm = rep(NA, nrow(rco@output))
	rcount.prof.segm[which.ok] = res
	rco@output$segmented = rcount.prof.segm
	rco
}

applyThresholds = function(rco, thresholds, min.CN) {
	if (length(rco@tumour) == 0) stop("Run readSeqsumm first!")
	if (length(rco@y) == 0) stop("Run applyFilters first!")
	if (ncol(rco@output) == 0) stop("Run runSeqnorm first!")
	if (ncol(rco@output) < 4) stop("Run runGLAD first!")
	cu = cut(rco@output$segmented, breaks=c(-Inf,thresholds,Inf), labels=min.CN+c(0:length(thresholds)))
	rco@output$CN = cu
	rco@thr = thresholds
	rco@minCN = min.CN
	rco
}

writeCNProfile = function(rco, folder) {
	if (ncol(rco@output) == 0) stop("Run runSeqnorm first!")
	write.table(rco@output, file=paste(folder,"/seqCNA_out.txt",sep=""), quote=FALSE,row.names=FALSE,sep="\t")
}

plotCNProfile = function(rco, folder=NULL) {
	if (ncol(rco@output) == 0) stop("Run runSeqnorm first!")
	profiles = rco@output
	if (.folderExists(folder)) {
		jpeg(file.path(folder,"seqCNA_out.jpg"),1920,1080)
		cex = 2.5
	} else {
		cex = 1.4
	}
	par(mfrow=c(1,1))
	par(mar=c(5,4,4,0.5))
	par(oma=c(0,0,0,0))
	binsX = 500
	bandX = ceiling(length(profiles$normalized) / binsX)
	binsY = 200
	bandY = diff(quantile(profiles$normalized[is.finite(profiles$normalized)],c(0.01,0.99))) / binsY
	suppressMessages(suppressWarnings(
	smoothScatter(postPlotHook=NULL, colramp=colorRampPalette(c("white","black")),
		bandwidth=c(bandX,bandY), nbin=c(binsX,binsY),
		profiles$normalized, cex=0.5, pch=19, col="#00000055",
		ylim=quantile(profiles$normalized[is.finite(profiles$normalized)],probs=c(0.005,0.995)),
		xlab="Genomic window index", ylab="Normalized ratio", frame.plot=FALSE,cex.axis=1.5,cex.lab=1.5)
	))
	abline(v=which(diff(as.numeric(as.factor(rco@seq)))!=0), lwd=2)
	minh = floor(min(profiles$normalized,na.rm=T))
	maxh = ceiling(max(profiles$normalized,na.rm=T))
	rang = maxh-minh
	if (rang>5) byh = 1
	else if (rang>1) byh = 0.2
	else if (rang>0.2) byh = 0.05
	else byh = 0.01
	abline(h=seq(minh,maxh,by=byh),lty=2,col="grey",lwd=2)
	if ("segmented"%in%colnames(profiles)) {		
		wh.ok = which(!is.na(profiles$segmented))
		wh.ch = wh.ok[diff(profiles$segmented[wh.ok])!=0]
		en = c(wh.ch, tail(wh.ok,1))
		st = c(1, en[1:(length(en)-1)]+1)
		vals = profiles$segmented[en]
		dummy = sapply(1:length(en), function(o) lines(c(st[o],en[o]),rep(vals[o],2),lwd=4,col="green"))
	}
	if ("CN"%in%colnames(profiles)) {
		for (cn in 0:length(rco@thr)) {
			ypos = median(profiles$segmented[profiles$CN==rco@minCN+cn], na.rm=T)
			text(0, ypos, paste("CN",rco@minCN+cn,sep=""), cex=cex, col="#ff0000", adj=1)
		}
		for (i in 1:length(rco@thr))
			abline(h=rco@thr[i], col="#ff0000", lwd=3)
	}
	if (.folderExists(folder))
		dev.off()
}

##################################################################
### SEQNORM FUNCTION
### this is the normalization function, which admits GC-based and
### paired normal modes, with different regression strategies
##################################################################
.seqnorm = function(x,y, nproc=1, resolution=1, use.normal=FALSE, skip=NULL, chr=rep(1,length(x)), 
		method=c("loess","cubic","quadratic")[3], out.dir=getwd(), out.jpg, plots=TRUE,
		lambdabreak=8, minRegions=7, maxRegions=35) {
	if (length(x) != length(y)) {
		if (! use.normal)
			stop("Tumour sample and GC content lengths differ.")
		else
			stop("Tumour and normal sample lengths differ.")
	}
	if (any(x[is.finite(x)]<0)) {
		if (! use.normal)
			stop("It seems that some GC content values are negative. Try recoding from -1 to NA in the GC file.")
		else
			stop("Hmmm... Something is not right, some read counts are negative!")
	}
	if (any(y[is.finite(y)]<0))
		stop("Hmmm... Something is not right, some read counts are negative!")
	if (nproc < 1)
		stop("Try setting nproc to a positive number. At least one processor core is needed!")
	if (resolution <= 0 || resolution >1)
		stop("resolution should be greater than 0 and as big as 1.")
	if (nproc > 64)
		warning("Warning: unless you really have a computer with more than 64 processor cores, try setting nproc to a lower number.", immediate.=TRUE)

	# work only with positive read counts and restore the rest after processing
	ret = rep(NA, length(x))
	xy.a = intersect(which(!is.na(x)), which(!is.na(y)))
	xy.a = intersect(xy.a, intersect(which(x>0),which(y>0)))
	if (! is.null(skip))
		xy.a = setdiff(xy.a, skip) # filter certain windows expected to misbehave
	x.a = x
	y.a = y
	chr.a = chr
	x = x[xy.a]
	y = y[xy.a]
	chr = chr[xy.a]
	chrtab = rle(as.vector(chr))$values
	chrlist = lapply(chrtab, function(chri) which(chr==chri))
	chrlist = chrlist[order(sapply(chrlist,head,1))]
	
	# work with bigger summarized windows if resolution < 1
	summ = 1
	if (resolution < 1) {
		chrlist.full = chrlist
		chrlist = list()
		max.ind = 0
		x.full = x
		y.full = y
		summ = floor(1/resolution)
		x = y = c()
		for (l in chrlist.full) {
			xchr = x.full[l]
			ychr = y.full[l]
			len = ceiling(length(xchr)/summ)
			x = c(x, sapply(split(xchr, rep(c(1:len),rep(summ,len))[1:length(xchr)]), mean, na.rm=T))
			y = c(y, sapply(split(ychr, rep(c(1:len),rep(summ,len))[1:length(ychr)]), mean, na.rm=T))
			chrlist[[length(chrlist)+1]] = c(1:len) + max.ind
			max.ind = tail(chrlist[[length(chrlist)]], 1)
		}
	}

	message("Creating segments...")
	flush.console()
	mins = sapply(chrlist,head,1)
	maxs = sapply(chrlist,tail,1)
	reg.method = ifelse(method=="loess", "loess", "lm")
	classicLomo = eval(call(reg.method, .makeFormula("x","y",method,use.normal,environment())))
	classic.rcount.prof = classicLomo$residuals/classicLomo$fitted
	cluster = makeCluster(rep("localhost",nproc), type="SOCK") 
	registerDoSNOW(cluster)
	i = 0
	bps = foreach(i=1:length(chrlist), .packages="GLAD") %dopar% {
		prof = classic.rcount.prof[chrlist[[i]]]
		input = cbind(c(1:length(prof)), prof, rep(1,length(prof)))
		colnames(input) = c("PosOrder","LogRatio","Chromosome")
		profileCGH = as.profileCGH(data.frame(input))
		res = glad(profileCGH, verbose=FALSE, bandwidth=1)$profileValues$Smoothing
		which(diff(res)!=0) + mins[i]-1
	}
	breaks = round(sort(c(1, unlist(bps), maxs[-length(maxs)], length(y)+1)))
	segmLength = diff(breaks)
	message("  Segments total: ", appendLF=FALSE)
	message(length(segmLength))
	flush.console()
	segmDomain = sapply(c(1:(length(breaks)-1)), function(i) { 
			Reduce("-",quantile(x[breaks[i]:(breaks[i+1]-1)],c(0.95,0.05)))
		} )
	segmRange = sapply(c(1:(length(breaks)-1)), function(i) { 
			Reduce("-",quantile(classic.rcount.prof[breaks[i]:(breaks[i+1]-1)],c(0.95,0.05)))
		} )
	# filter regions based on three different metrics
	minDomain = 0.95 # minimum span over X of a segment's middle 90% quantile (in quantiles)
	maxRange = 0.85 # maximum span over Y of a segment's middle 90% quantile (in quantiles)
	minLength = 0.85 # minimum segment length
	stepp = 0.025
	minRegions = 7
	maxRegions = 35
	.getSegmLm = function() {
		segmLm = intersect(which(segmLength>quantile(segmLength,minLength)), which(segmRange<quantile(segmRange,maxRange)))
		intersect(segmLm, which(segmDomain>quantile(segmDomain,minDomain)))
	}
	segmLm = .getSegmLm()
	while (length(segmLm) < min(minRegions, length(breaks)-1) && maxRange < 1-stepp) {
		minDomain = minDomain - stepp
		maxRange = maxRange + stepp
		minLength = minLength - stepp
		segmLm = .getSegmLm()
	}
	while (length(segmLm) > maxRegions && minDomain < 1-stepp) {
		minDomain = minDomain + stepp
		maxRange = maxRange - stepp
		minLength = minLength + stepp
		segmLm = .getSegmLm()
	}
	message("  Segments passing QC: ", appendLF=FALSE)
	message(length(segmLm))
	flush.console()
	classic.rcount.prof = scale(classicLomo$residuals/classicLomo$fitted)
	if (length(segmLm) < min(minRegions, length(breaks)-1)) {
		message("  Not enough segments, falling back to typical regression.")
		flush.console()
		ret[xy.a] = classic.rcount.prof
		return(ret)
	}
	
	message("Fitting regression models...")
	flush.console()
	if (resolution < 1) sel = seq(1,length(x.full),l=100000)
	xminwh = which.min(x)
	xmin = x[xminwh]
	ymin = y[xminwh]
	xmaxwh = which.max(x)
	xmax = x[xmaxwh]
	ymax = y[xmaxwh]
	fit = foreach(i=c(1:(length(breaks)-1))[segmLm], .export=c(".makeFormula")) %dopar% {
		rang = breaks[i]:(breaks[i+1]-1)
		xLm = x[rang]
		yLm = y[rang]
		if (reg.method=="loess") {
			xLm = c(xLm, xmin,xmax)
			yLm = c(yLm, ymin,ymax)
		}
		lomo = eval(call(reg.method, .makeFormula("xLm","yLm",method,use.normal,environment())))
		if (resolution < 1) p = x.full[sel]
		else p = x
		predict(lomo, data.frame("xLm"=p))
	}
	fit = do.call(rbind, fit)
	# subtract shift to ref regression
	baseLomo = median(classicLomo$fitted)
	fit = t(apply(fit, 1, function(i) i / (median(i)/baseLomo))) 
	# for each point, merge sub-regression estimations
	chunk.len = ceiling(ncol(fit)/nproc)
	n = 0
	medianFit = foreach(n=1:nproc) %dopar% {
		fits = c()
		for (i in c(1:chunk.len)+(n-1)*chunk.len)
			if (i <= ncol(fit))
				fits = c(fits, median(fit[,i]))
		fits
	}
	stopCluster(cluster)
	medianFit = unlist(medianFit)
	if (resolution < 1) {
		xSel = x.full[sel]
		lomo = eval(call(reg.method, .makeFormula("xSel","medianFit",method,use.normal,environment())))
		lomoFitted = predict(lomo, data.frame(xSel=x.full))
		rcount.prof = y.full / lomoFitted
	} else {
		lomo = eval(call(reg.method, .makeFormula("x","medianFit",method,use.normal,environment())))
		lomoFitted = lomo$fitted
		rcount.prof = y / lomoFitted
	}
	rcount.prof = scale(rcount.prof)

	if (plots) {
		message("Plotting normalization...")
		flush.console()
		if (.folderExists(out.dir)) {
			jpeg(file.path(out.dir,"seqCNA_segment_QC.jpg"),1000,700)
			segmLengthCut = cut(segmLength, breaks=c(-Inf,quantile(segmLength,0.5),Inf), labels=c(1,2))
			layout(matrix(c(1:4),ncol=2), widths=c(1,4), heights=c(4,1))
			par(oma=c(0,0,0,0))
			par(mar=c(4,0,1,1))
			d = density(segmRange[which(segmLengthCut==2)], from=min(segmRange), to=max(segmRange), n=10000)
			plot(-d$y, d$x, frame.plot=FALSE, xaxt="n",yaxt="n", type="p", col="black", xlab="",ylab="", pch=19, cex=0.25)
			d = density(segmRange[which(segmLengthCut==1)], from=min(segmRange), to=max(segmRange), n=10000)
			points(-d$y, d$x, col="grey", pch=19, cex=0.25)
			par(mar=c(1,0,1,1))
			plot.new()
			legend("topleft",c("Short segments","","Long segments"), col=c("Grey","White","Black"),lwd=4, box.col="white")
			par(mar=c(4,4,1,1))
			plot(segmDomain, segmRange, col=c("#88888888","#000000aa")[segmLengthCut],
				frame.plot=FALSE, xlab="Segment domain", ylab="Segment range", cex=1, pch=19)
			abline(h=quantile(segmRange,0.75))
			abline(v=quantile(segmDomain,0.75))
			par(mar=c(0,4,1,1))
			d = density(segmDomain[which(segmLengthCut==1)], from=min(segmDomain), to=max(segmDomain), n=10000)
			plot(d$x, d$y, frame.plot=FALSE, xaxt="n",yaxt="n", type="p", xlab="",ylab="", col="grey", pch=19, cex=0.25)
			d = density(segmDomain[which(segmLengthCut==2)], from=min(segmDomain), to=max(segmDomain), n=10000)
			points(d, col="black", pch=19, cex=0.25)
			dev.off()
		}
		if (.folderExists(out.dir)) {
			jpeg(out.jpg,2000,1600)
			cex = 2
		} else {
			cex = 1
		}
		par(mfrow=c(1,2))
		par(oma=c(0,0,0,0))
		par(mar=c(4,5,1,1))
		qu = quantile(rcount.prof, c(0.001,0.999))
		wh = which(rcount.prof>qu[1] & rcount.prof<qu[2])
		plot(density(rcount.prof[wh],n=5000,adjust=1/summ), col="#ff6600ff", xlim=quantile(rcount.prof, c(0.01,0.99)), main="",
			xlab="Normalized ratio", frame.plot=FALSE,cex.axis=cex,cex.lab=cex, lwd=2)
		if (resolution == 1) {
			lines(density(classic.rcount.prof[wh],n=5000), col="#0000FFFF", lwd=2)
		}
		xlab = ifelse(use.normal, "Window counts (normal)", "GC content")
		suppressMessages(suppressWarnings(
		smoothScatter(colramp=colorRampPalette(c("white","black")),
			x,y, col="#00000022", cex=1, pch=19, xlab=xlab, ylab="Window counts (tumour)",
			xlim=quantile(x, c(0.001,0.999)), ylim=quantile(y, c(0.01,0.99)), cex.axis=cex, cex.lab=cex)
		))
		if (resolution < 1) {
			ord = order(x.full)[sel]
			lines(x.full[ord], lomoFitted[ord], pch=19, col="#FF6600FF", lwd=3)
		} else {
			ord = order(x)
			lines(x[ord], lomoFitted[ord], pch=19, col="#FF6600FF", lwd=3)
			lines(x[ord], classicLomo$fitted[ord], pch=19, col="#0000FFFF", lwd=3)
			legend("topright",c("Standard regression","","seqnorm regression"),
				col=c("#0000FFFF","#FFFFFF00","#FF6600FF"),lwd=4, box.col="#00000000", cex=cex*1.2)
		}
		if (.folderExists(out.dir))  dev.off()
	}

	ret[xy.a] = rcount.prof
	ret
}

##################################################################
### AUXILIARY FUNCTIONS
##################################################################
.folderExists = function(x) {
	if (is.null(x)) return(0)
	if (is.na(x)) return(0)
	x = gsub("[\\/]*$","",x)
	if (!file.exists(x)) return(0)
	ret = 1
	wd = getwd()
	try.wd = try(setwd(wd), silent=TRUE)
	if (class(try.wd) == "try-error")
		ret = 0
	else
		setwd(wd)
	ret
}
.outputExists = function(x) length(list.files(path=x, pattern="^seqsumm_out.txt$")) > 0
.resample = function(info, new.win, nproc=2) {
	new.win = new.win*1000
	win = diff(info$win.start[1:2])
	if (new.win %% win != 0 || new.win < win)
		stop(paste("New window size should be an exact multiple of ",win/1000,".",sep=""))
	m.factor = new.win / win
	info$reads.gc[info$reads.gc==0] = NA
	info$reads.mapq[info$reads.mapq==0] = NA
	chrs = rle(as.character(info$chrom))$values
	cluster = makeCluster(rep("localhost",nproc), type="SOCK") 
	registerDoSNOW(cluster)
	chr = 0
	ss = foreach(chr=chrs) %dopar% {
		s = list()
		chr.info = info[info$chrom==chr, ]
		l = ceiling(nrow(chr.info)/m.factor)
		spls = lapply(c(3:ncol(info)), function(i) split(chr.info[,i], rep(c(1:l),rep(m.factor,l))[1:nrow(chr.info)]))
		for (i in 1:2) # weighted average for GC content and mapping quality
			s[[i]] = sapply(c(1:length(spls[[i]])), function(j) 
				sum(spls[[i]][[j]] * spls[[3]][[j]], na.rm=TRUE) /  sum(spls[[3]][[j]], na.rm=TRUE))
		for (i in 3:(ncol(info)-2)) # sum for reads
			s[[i]] = sapply(c(1:length(spls[[i]])), function(j) 
				sum(spls[[i]][[j]], na.rm=TRUE))
		do.call(cbind, s)
	}
	stopCluster(cluster)
	dat = do.call(rbind, ss)
	l = sapply(ss, nrow)
	ch = rep(chrs, l)
	win.st = unlist(lapply(l, function(i) c(0:(i-1))*new.win))
	newinfo = data.frame(chrom=ch, win.start=win.st, reads.gc=dat[,1], reads.mapq=dat[,2], stringsAsFactors=FALSE)
	for (i in 3:ncol(dat))
		newinfo = cbind(newinfo, dat[,i])
	colnames(newinfo) = colnames(info)
	newinfo
}
.makeFormula = function(xn, yn, method, use.normal, env) {
	formulae = paste(yn,"~",xn)
	if (method != "loess") formulae = paste(formulae, "+I(",xn,"^2)+I(",xn,"^3)")
	if (method == "quadratic") formulae = paste(formulae, "+I(",xn,"^4)")
	if (use.normal) formulae = paste(formulae, "-1")
	as.formula(formulae, env=env)
}
.buildGenomeInfo = function(info, info_len, chrs, win, nproc=2) {
	csum = cumsum(ceiling(info_len$length/1000))
	names(csum) = info_len$chr
	en = csum[chrs]
	st = c(1, en[1:(length(en)-1)]+1)
	names(st) = names(en)
	cluster = makeCluster(rep("localhost",nproc), type="SOCK") 
	registerDoSNOW(cluster)
	chr = 0
	ss = foreach(chr=chrs) %dopar% {
		s = list()
		chr.info = info[st[chr]:en[chr], ]
		l = ceiling(nrow(chr.info)/win)
		for (i in 1:3) {
			if (win == 1) {
				s[[i]] = chr.info[,i]
			} else {
				spl = split(chr.info[,i], rep(c(1:l),rep(win,l))[1:nrow(chr.info)])
				s[[i]] = sapply(spl, mean, na.rm=T)
			}
		}
		do.call(cbind, s)
	}
	stopCluster(cluster)
	dat = do.call(rbind, ss)
	l = sapply(ss, nrow)
	ch = rep(chrs, l)
	win.st = unlist(lapply(l, function(i) c(0:(i-1))*win))
	newinfo = cbind(ch, win.st, dat)
	colnames(newinfo) = c("Chr", "Start", colnames(info))
	newinfo = as.data.frame(newinfo)
	for (i in 3:5)
		newinfo[,i] = as.numeric(as.character(newinfo[,i]))
	newinfo
}
.autoTrim = function(rc, from=0.98, by=0.0005, nproc=2) {
	to = 1 - 10 / length(rc)
	quants = seq(from,to,by)
	cluster = makeCluster(rep("localhost",nproc), type="SOCK") 
	registerDoSNOW(cluster)
	chunk.len = ceiling(length(quants)/nproc)
	n = 0
	ws = foreach(n=1:nproc, .packages="adehabitatLT") %dopar% {
		w = c()
		for (i in c(1:chunk.len)+(n-1)*chunk.len)
			if (i <= length(quants))
				w = c(w, wawotest(rc > quantile(rc, quants[i]))[1])
		w
	}
	stopCluster(cluster)
	ws = unlist(ws)
	names(ws) = quants
	ws.smoothed = loess(ws~quants, span=0.5)$fitted
	whm = which.max(ws.smoothed)
	rang = whm:length(ws.smoothed)
	ws.smoothed = ws.smoothed[rang]
	y = ws.smoothed - min(ws.smoothed)
	y = y / max(y)
	x = quantile(rc,quants)[rang]
	x = x - min(x)
	x = x / max(x)
	last.valley = length(y) - tail(rle(diff(y)<0)$lengths,1)
	sel = last.valley : length(y)
	if (all(diff(y[sel])>=0))
		sep = length(quants)
	else
		sep = which.min(x[sel]+y[sel]) + whm-1 + last.valley-1
	trim = quants[sep]
	trim
}
