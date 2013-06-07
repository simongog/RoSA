source("../src/basic_functions.R")
library("tikzDevice")




fn="web-4GB"
data <- data_frame_from_key_value_pairs( paste("query_performance_details_",fn,".txt", sep="" ) )

data['rtime_full'] <- data['rtime_full']/data['full_queries']
data['utime_full'] <- data['utime_full']/data['full_queries']

data <- subset(data, data[["pattern_file_name"]]=="../pattern/web-4GB.40.1000.0.75.125.pattern")
#data <- subset(data, data[["pattern_file_name"]]=="../pattern/web-4GB.40.1000.0.8.12.pattern")
data <- subset(data, data[["fac_dens"]] < 513)

d_mean <- aggregate(data[c('rtime_full','utime_full','fac_dens','disk_access_per_query','phase')] ,list(data[['fac_dens']],data[['phase']]),mean)

dd <- d_mean[order(d_mean['fac_dens'],d_mean['phase']),]
phases=6
tab <- dd[['rtime_full']]
dim(tab) <- c(phases, length(tab)/phases)
for(i in seq(phases, 2)){
	tab[i,] <- tab[i,]-tab[i-1,]
}

#mycol <- topo.colors(2*phases)[seq(1,2*phases,2)]
mycol <- terrain.colors(phases)

fac_denss <- unique(dd[['fac_dens']])
fac_dens_label <- paste(rep("$",length(fac_denss)), fac_denss,rep("$",length(fac_denss)),sep="")

tikz("/Users/sgog/Downloads/cpm/fig/detailed_query_time.tex", width="2.4", height="2.3")

par( oma=c(0.1,0.1,0.1,0.1) )
par( mar=c(1.5,1.5,0.1,0.1) )

barplot(tab, names.arg=rep("", length(fac_dens_label)), col=mycol, axes=F)

legend("topleft", legend=rev(c("internal matching","load disk block","build block tree","block tree matching","load text","match pattern")),
		fill=rev(mycol), bty="o", box.lwd=0, bg="white", inset=c(0.01,-0.05), y.intersp=0.8 )

#grid(col="gray")

barp <- barplot(tab, ylab="Runtime per query phase",col=mycol, add=T, cex.axis=0.8, yaxt="n")
axis(1, at=barp, seq(1, length(fac_dens_label)), labels=fac_dens_label, cex.axis=0.9, lty=0, line=-0.6)
axis(2, line=-0.5, cex.axis=0.8, lty=0)
axis(2, line=-0.2, cex.axis=0.8, labels=F )

dev.off()



tikz("/Users/sgog/Downloads/cpm/fig/detailed_space.tex", width="2.4", height="2.3")

par( oma=c(0.1,0.1,0.1,0.1) )
par( mar=c(1.5,1.5,0.1,0.1) )



data <-  data_frame_from_key_value_pairs("../space_usage/space_web-4GB.txt")
data <- data[order(data[['fac_dens']]),]

data <- subset(data, data[['fac_dens']] %in% c(0,1,4,16,64,256))

size <- c( data[['compact_text_in_megabyte']], data[['header_in_megabyte']], data[['bp_ct_in_megabyte']], data[['lcp_in_megabyte']], data[['sa_in_megabyte']] )

dim(size) <- c( nrow(data), 5)
size <-	t(size)

rcol=terrain.colors(5)

barp <- barplot( size/(data[["n"]]/(1024*1024)), col=rcol,
	     ylab="external index size in percent of orig. text", ylim=c(0, 3.2), yaxt="n"   )

axis(1, at=barp, labels = data[['fac_dens']], cex.axis=0.9, lty=0, line=-0.6)
axis(2, line=-0.5, cex.axis=0.8, lty=0)
axis(2, line=-0.2, cex.axis=0.8, labels=F)

legend("topright", legend=rev(c("text representation","header","topology","LCP","text pointer")), fill=rev(rcol),
	   inset=c(-0.1,-0.05),
	   box.lwd=0, y.intersp=0.8)

dev.off()
