source("../src/basic_functions.R")

data <-  data_frame_from_key_value_pairs("space_web-4GB.txt")

data <- data[order(data[['fac_dens']]),]


size <- c( data[['compact_text_in_megabyte']], data[['header_in_megabyte']], data[['bp_ct_in_megabyte']], data[['lcp_in_megabyte']], data[['sa_in_megabyte']] )

dim(size) <- c( nrow(data), 5)
size <-	t(size)

rcol=terrain.colors(6)

barplot( size/(data[["n"]]/(1024*1024)), names.arg=data[['fac_dens']], col=rcol,
	     ylab="external index size in percent of orig. text", xlab="R", main="space breakdown web-4GB"   )

legend("topright", legend=c("text representation","header","topology","LCP","text pointer"), fill=rcol)
