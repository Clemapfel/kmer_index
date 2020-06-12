# Created by: clem
# Created on: 6/12/20

filter_extremes <- function(data, interval=0.975) {

  # hash for
  data_for <- data[data$name == "hash_for",]
  mean_for <- mean(data_for$real_time)
  stddev_for <- sd(data_for$real_time)
  error_for <- qnorm(interval)*stddev_for/sqrt(nrow(data_for))
  lower_for <- mean_for - error_for
  upper_for <- mean_for + error_for

  # hash fold
  data_fold <- data[data$name == "hash_fold",]
  mean_fold <- mean(data_fold$real_time)
  stddev_fold <- sd(data_fold$real_time)
  error_fold <- qnorm(interval)*stddev_fold/sqrt(nrow(data_fold))
  lower_fold <- mean_fold - error_fold
  upper_fold <- mean_fold + error_fold

  return(data[(data$name == "hash_for" & data$real_time < upper_for & data$real_time > lower_for)
                |(data$name == "hash_fold" & data$real_time < upper_fold & data$real_time > lower_fold),])
}

library(ggplot2)

setwd("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/map_vs_map/")
data <- read.csv("2020-06-12_18-36-59.csv")

filter = 0.999;

# filter out gbenchmark averages
data_at <- data[data$name == "da_at" | data$name == "std_at" | data$name == "robin_hood_at",]

plot_at <- qplot(name, real_time, data=data_at, fill=name, geom="boxplot")
plot_at = plot_at + ggtitle(label="unordered_map::at(size_t) performance for different map implementations",
                            subtitle="TODO")
plot_at = plot_at + theme(plot.title=element_text(face="bold"))

data_insert <- data[data$name == "da_insert" | data$name == "std_insert" | data$name == "robin_hood_insert",]

plot_insert <- qplot(name, real_time, data=data_insert, fill=name, geom="boxplot")
plot_insert = plot_insert + ggtitle(label="unordered_map::insert(size_t, uint32_t) performance for different map implementations",
                                    subtitle="TODO")
plot_insert = plot_insert + theme(plot.title=element_text(face="bold"))
