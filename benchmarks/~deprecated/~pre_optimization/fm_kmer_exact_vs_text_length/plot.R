# plot search k vs query size
setwd("~/Documents/Workspace/kmer_index/source/benchmarks/fm_kmer_exact_vs_text_length")
data_raw = read.csv("2020-05-07_17-26-48.csv")
data_kmer = data_raw[data_raw$name == "multi_10_11_12_kmer_search_mean",]
data_fm = data_raw[data_raw$name == "fm_search_mean",]

par(mfrow=c(2,1))

plot(data_kmer$text_size, data_kmer$real_time, main="exact search", xlab="text size", ylab="runtime (ns)", type="l", col="red");
plot(data_fm$text_size, data_fm$real_time, col="green", type = "l")

