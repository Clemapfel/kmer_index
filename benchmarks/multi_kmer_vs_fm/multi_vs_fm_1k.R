# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 6/19/20

library(ggplot2)

setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/")
data <- read.csv("final_1k_3k_2020-06-24T00-13-32+02-00.csv")
data_fm = data[data$name == "fm_median",]

#data <- read.csv("TODO:")
data_multi = data[data$name == "multi_kmer_median",]

plot_diff <- ggplot()
diff = data_multi$real_time - data_fm$real_time
diff_log = sign(diff) * log(diff)
plot_diff = plot_diff + geom_segment(mapping=aes(x=data_multi$query_length, y=0, xend=data_multi$query_length, yend=diff_log))
#plot_diff = plot_diff + #coord_cartesian(ylim=c(-5000, 10000))
plot_diff
