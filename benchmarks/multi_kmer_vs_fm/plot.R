# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 6/17/20

library(ggplot2)

setwd("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/")
data <- read.csv("2020-06-16_19-59-09.csv")

get_geom_line <- function(name, color) {
  return(geom_line(mapping=aes(query_length, real_time, color=name), data=data[data$name==name,]))
}

plot = ggplot() + get_geom_line("multi_kmer_mean") + get_geom_line("single_kmer_mean") + get_geom_line("fm_median")
plot = plot + ggtitle("search performance", subtitle="text size = 1e6, queries randomized each call, 3 benchmark cycles per query size")
plot = plot + theme(plot.title=element_text(face="bold"))

labels = c("fm_index", "kmer_index<20, TODO>", "kmer_index<20>")
colors = c("red", "darkorchid3", "blue3")
plot = plot + scale_color_manual(name="index", labels=labels, values=colors)

# add horizontal lines for each k
k <- 20
for (i in 1:10) {

  plot = plot + geom_vline(xintercept=k*i, color="darkgrey", linetype="dotted", size=1)
}

plot = plot + scale_x_continuous(breaks=seq(0, 200, 10))




plot