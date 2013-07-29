setClass("SeqCNAInfo", representation(tumour="list", normal="list",
							seq="character", pos="numeric", build="character", win="numeric",
							x="numeric", y="numeric",
							skip="numeric", output="data.frame",
							thr="numeric", minCN="numeric",
							gc="numeric", mapq="numeric"))