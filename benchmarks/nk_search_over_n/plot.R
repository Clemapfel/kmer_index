# plot search k vs query size
#setwd("~/Documents/Workspace/kmer_index/source/benchmarks/nk_search_over_k/")
data = read.csv("2020-05-22_23-02-10.csv")

k <- 6

# filter sdev
data <- data[data$name=="kmer_search_mean",]

plot(data$query_length, data$real_time, type="l")
for (i in 1:10)
  abline(v=i*k, col="red")




