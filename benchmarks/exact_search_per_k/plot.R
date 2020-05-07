# plot search k vs query size
setwd("~/Documents/Workspace/kmer_index/source/benchmarks/exact_search_per_k")
data = read.csv("2020-05-06_16-29-02.csv")

min_k = 100;
max_k = 0;

# filter sdev
data = data[data$name=="search_mean",]

different_k = unique(data$k)
par(mfrow=c(1, length(different_k)))

for (k in different_k) {

  current = data[data$k == k,]
  plot(current$query_size, current$real_time, main=paste("k = ", k), type="l")

  for (i in 1:4)
    abline(v=i*k, col = "red")

}


