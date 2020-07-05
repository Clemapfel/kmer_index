# Created by: clem
# Created on: 26.06.20

library(ggplot2)
library(ggpubr)
library(dplyr)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/just_k/")
data = read.csv("lncrna_complete.csv")
data = data[!grepl("stddev", data$name, fixed=TRUE)  & data$query_length >=5,]

fm_color_label = "fm"
kmer_color_label = "kmer"

smooth_size = 3

get_plot = function(text_length)  {


  data_cur = data[data$text_length == text_length,]
  fm = data_cur[data_cur$name=="fm_median",];  #sic
  kmer = data_cur[data_cur$name=="kmer_median",]

  plot = ggplot() + geom_line(mapping=aes(x=query_length, y=real_time, color=name), data=fm)
  plot = plot + geom_smooth(mapping=aes(x=query_length, y=real_time, color=name), method="lm", data=fm, size=2)

  plot = plot + geom_line(mapping=aes(x=query_length, y=real_time, color=name), data=kmer)
  plot = plot + geom_smooth(mapping=aes(x=query_length, y=real_time, color=name), method="lm", data=kmer, size=2,)

  plot = plot + coord_cartesian(ylim=c(min(kmer$real_time), max(fm$real_time)))
  plot = plot + scale_x_continuous(breaks=seq(1, 30, 1))
  plot = plot + ggtitle(label=paste("text length =", text_length))

  return(plot)
}

get_plot(1e6)










ggsave("kmer_search_fm_vs_single.png", plot_e6)