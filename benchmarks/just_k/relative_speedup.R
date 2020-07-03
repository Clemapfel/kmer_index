library(ggplot2)


# to 1e6
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/just_k/")
data_big = read.csv("to_1e7_complere_2020-07-02T18-27-15+02-00.csv")
data_big = data_big[!grepl("stddev", data_big$name, fixed=TRUE) ,]

for (text_length in sort(unique(data_big$text_length))) {

  data_cur = data_big[data_big$text_length == text_length,]
  fm = data_cur[data_cur$name=="fm_median",]  #sic
  kmer = data_cur[data_cur$name=="kmer_median",]
  print(paste(text_length, mean(fm$real_time / (fm$real_time - kmer$real_time))))
}

#####################

print_plot = function(text_length) {

    fm = data_big[data_big$name == "fm_median" & data_big$text_length == text_length,]
    kmer = data_big[data_big$name == "kmer_median" & data_big$text_length == text_length,]

    plot = ggplot()
    plot = plot + geom_segment(aes(x=kmer$query_length, xend=kmer$query_length, y=0, yend=kmer$real_time - fm$real_time))
    plot = plot + scale_x_continuous(breaks=seq(1, 30, 1)) + ggtitle(label="kmer - fm", subtitle=paste("text size = ", text_length))
    print(plot)
}
