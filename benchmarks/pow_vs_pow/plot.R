# Created by: clem
# Created on: 6/8/20

library(ggplot2)

setwd("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/pow_vs_pow/")
data <- read.csv("2020-06-09_16-08-44.csv")
data_filtered <- data[data$name=="std_pow" | data$name=="fast_pow" | data$name=="recursive_pow" | data$name=="bit_pow",] # filter gbenchmark averages
#data_std <- data[data$name=="std_pow",]
#data_fast <- data[data$name=="fast_pow",]

plot <- qplot(name,
              real_time,
              data = data_filtered,
              geom = "boxplot",
              notch = TRUE,
              xlab = "",
              ylab = "runtime (ns)",
              ylim = c(110, 175),
              main = "runtime for different pow implementations",
              fill = name)

plot
