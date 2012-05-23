

# Read a file called file_name and create a data frame the following way
# (1) Parse all the lines of the form
# '# key = value'
# (2) Each unique key gets a column
data_frame_from_key_value_pairs <- function(file_name){
	lines <- readLines(file_name)
	lines <- lines[grep("^#.*=.*",lines)]
	d <- gsub("^#","",gsub("[[:space:]]","",unlist(strsplit(lines,split="="))))
	keys <- unique(d[seq(1,length(d),2)])
	keynr <- length(keys)
	dd <- d[seq(2,length(d),2)]
	dim(dd) <- c( keynr, length(dd)/keynr )
	data <- data.frame(t(dd))	
	names(data) <- keys	

	for (col in keys){
		t <- as.character(data[[col]])
		suppressWarnings( tt <- as.numeric(t) )
		if ( length( tt[is.na(tt)] ) == 0 ){ # if there are not NA in tt
			data[[col]] <- tt 
		}
	} 	
	data
}

format_str_fixed_width <- function(x, width=4){
	sx <- as.character(x)
	if ( nchar(sx) < width ){
		for (i in 1:(width-nchar(sx))){
			sx <- paste("\\D",sx, sep="")
		}
	}
	sx
}
