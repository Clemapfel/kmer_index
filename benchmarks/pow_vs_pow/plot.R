# Created by: clem
# Created on: 6/8/20

# plot runtime for pows
library(ggplot2)

setwd("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/pow_vs_pow/")
data <- read.csv("2020-06-09_16-52-16.csv")

# filter out gbenchmark averages
data_filtered <- data[data$name=="switch_pow"
                        | data$name=="x*x"
                        | data$name=="recursive_pow"
                        | data$name=="bit_pow",]

plot <- qplot(name,
              log(real_time),
              data = data_filtered,
              geom = "boxplot",
              notch = TRUE,
              xlab = "",
              ylab = "log(runtime)",
              #ylim = c(90, 230),
              main = "runtime distribution for different pow implementations (base and exponent randomized each call)",
              fill = name)

ggsave("pow_vs_pow.png", plot,  width=30, height=15, units="cm")
