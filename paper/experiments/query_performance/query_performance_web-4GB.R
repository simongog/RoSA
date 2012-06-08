library(xtable) # if not installed call install.packages("xtable")

source("../../../src/basic_functions.R")

data <- data_frame_from_key_value_pairs( "query_performance_web-4GB.txt" )

full_queries <- data[['full_queries']]
if ( max(full_queries) != min(full_queries) ){
	cat("WARNING: Thre are different values for 'full_queries' in the file\n")
}
full_queries <- full_queries[1]

in_memory_queries <- data[['in_memory_queries']]
if ( max(in_memory_queries) != min(in_memory_queries) ){
	cat("WARNING: Thre are different values for 'in_memroy_queries' in the file\n")
}
in_memory_queries <- full_queries[1]

# data <- split(data, data[["pattern_file_name"]])
disk_access <- aggregate( data['disk_access_per_query'], by = list( data[['pattern_file_name']]), FUN=mean )
rtime_full <- aggregate( data['rtime_full'], by = list( data[['pattern_file_name']]), FUN=mean )
utime_full <- aggregate( data['utime_full'], by = list( data[['pattern_file_name']]), FUN=mean )
rtime <- aggregate( data['rtime'], by = list( data[['pattern_file_name']]), FUN=mean )
avg_depth <- aggregate( data['avg_depth'], by = list( data[['pattern_file_name']]), FUN=mean )

disk_access_table <- c(rep(0,25))
dim(disk_access_table) <- c(5,5)

rtime_full_table <- c(rep(0,25))
dim(rtime_full_table) <- c(5,5)

utime_full_table <- c(rep(0,25))
dim(utime_full_table) <- c(5,5)

rtime_in_memory_table <- c(rep(0,25))
dim(rtime_in_memory_table) <- c(5,5)

avg_depth_table <- c(rep(0,25))
dim(avg_depth_table) <- c(5,5)

i <-  1	
lens <- c("4","10","20","40","100")
for ( len in lens ) {
	j <- 1
	for (range in c("1.2", "8.12","75.125","750.1250","7500.12500")){
		pattern <- paste(".*web-4GB\\.",len,"\\..*",range,".*\\.pattern",sep="")
		idx <- grep(pattern, disk_access[[1]]) 
		disk_access_table[i, j] <- disk_access[[2]][idx]
		rtime_full_table[i, j] <- rtime_full[[2]][idx]
		utime_full_table[i, j] <- utime_full[[2]][idx]
		rtime_in_memory_table[i, j] <- rtime[[2]][idx]
		avg_depth_table[i, j] <- avg_depth[[2]][idx]
		j <- j+1
	}
	i <- i+1
}

utime_full_table <- (utime_full_table / full_queries) * 1000
rtime_full_table <- (rtime_full_table / full_queries) * 1000
rtime_in_memory_table <- (rtime_in_memory_table / in_memory_queries) * 1000

print_latex <- function( f_table, f_digits, f_total_len = 0 ){
	latex_lens <- c("\\D\\D4","\\D10","\\D20","\\D40","100") 
	f_table <- round(f_table, f_digits)
	table <- apply( f_table, MARGIN=c(1,2), FUN=format_str_fixed_width, width=f_total_len )
	cr_table <- data.frame( table )
	cr_table <- cbind( rep("", nrow(cr_table)), cr_table)
	rownames( cr_table ) <- latex_lens
	print( xtable( cr_table, digits=f_digits ), type="latex", 
		   sanitize.rownames.function = identity, 
		   sanitize.text.function = identity 
		   )
}

print_latex( disk_access_table, 1 )
print_latex( rtime_full_table, 0, 4 )
print_latex( utime_full_table, 0, 4 )

