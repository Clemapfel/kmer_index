library(ggplot2)
library(tidyr)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/alphabet/")
data_raw = read.csv("alphabet_1e6_2020-07-03T18-22-20+02-00.csv")
data = data_raw[data_raw$name=="kmer_median",]

plot = ggplot()
x = unique(data$query_length)
plot = plot + geom_segment(mapping=aes(x=, xend=alphabet_size, y=0, yend=real_time, color=alphabet_size), data=data[data$query_length==5,])










plot