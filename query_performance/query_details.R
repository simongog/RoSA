source("../src/basic_functions.R")

fn="web-256MB"
data <- data_frame_from_key_value_pairs( paste("query_performance_details_",fn,".txt", sep="" ) )

data['rtime_full'] <- data['rtime_full']/data['full_queries']
data['utime_full'] <- data['utime_full']/data['full_queries']

d_mean <- aggregate(data[c('rtime_full','utime_full','fac_dens','disk_access_per_query','phase')] ,list(data[['fac_dens']],data[['phase']]),mean)

dd <- d_mean[order(d_mean['fac_dens'],d_mean['phase']),]

tab <- dd[['rtime_full']]
dim(tab) <- c(5, length(tab)/5)
for(i in seq(2,5)){
	tab[i] <- tab[i]-tab[i-1]
}

mycol=head(topo.colors(10),5)
barplot(tab, names.arg=unique(dd[['fac_dens']]), ylab="elapsed time per query phase in [millisec]",col=mycol, xlab="K",
		main=paste("Detailed runtime of a count query ", "l=40, k=100",sep=""))
legend("topleft", legend=rev(c("condensed BWT matching","load disk block","build block tree","load text","match pattern")),
		fill=rev(mycol))

