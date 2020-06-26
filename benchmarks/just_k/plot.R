# Created by: clem
# Created on: 26.06.20

library(ggplot2)
library(ggpubr)
library(dplyr)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/just_k/")
data_raw = read.csv("2020-06-25T23-47-49+02-00.csv")
data_raw = data_raw[!grepl("stddev", data_raw$name, fixed=TRUE) ,]

# difference per text_size
speedup = matrix(nrow=length(unique(data_raw$text_length)), ncol = 2)
colnames(speedup) = c("text_length", "speedup")

i = 0
for (text_length in unique(data_raw$text_length)) {

  speedup[i, "text_length"] = text_length
  data_cur = data_raw[data_raw$text_length == text_length,]
  fm = data_cur[data_cur$name=="single_kmer_median",]  #sic
  kmer = data_cur[data_cur$name=="kmer_single_median",]
  speedup[i, "speedup"] = mean(fm$real_time / (fm$real_time - kmer$real_time))
  i = i+1
}

speedup

title = "search runtime for kmer"
fm_offset = 0.1
fm_size = 5
kmer_size = 5

kmer_color_label = "kmer::kmer_index<k>"
fm_color_label = "seqan3::fm_index"





#1e3
########################################################################################
text_length = 1e3
data = data_raw[data_raw$text_length == text_length,]
data_kmer = data[data$name == "kmer_single_median",]
data_fm = data[data$name == "single_kmer_median",] #sic, typod the benchmark name

plot_e3 = ggplot()
plot_e3 = plot_e3 + geom_segment(data=data_fm, aes(x=query_length+fm_offset, xend=query_length+fm_offset, y=0, yend=real_time, color=fm_color_label), size=fm_size)
plot_e3 = plot_e3 + geom_segment(data=data_kmer, mapping=aes(x=query_length, xend=query_length, y=0, yend=real_time, color=kmer_color_label), size=kmer_size)

plot_e3 = plot_e3 + scale_x_continuous(name="k (query length)", breaks=seq(0, 30, 1)) + scale_y_continuous(name="runtime (ns)")
plot_e3 = plot_e3 + geom_hline(yintercept=0, color="grey")
plot_e3 = plot_e3 + ggtitle(label=title, subtitle=paste("text length =", text_length, "| mean of 150 benchmark cycles per k")) + theme(plot.title=element_text(face="bold"))

#1e4
########################################################################################
text_length = 1e4
data = data_raw[data_raw$text_length == text_length,]
data_kmer = data[data$name == "kmer_single_median",]
data_fm = data[data$name == "single_kmer_median",] #sic, typod the benchmark name

plot_e4 = ggplot()
plot_e4 = plot_e4 + geom_segment(data=data_fm, aes(x=query_length+fm_offset, xend=query_length+fm_offset, y=0, yend=real_time, color=fm_color_label), size=fm_size)
plot_e4 = plot_e4 + geom_segment(data=data_kmer, mapping=aes(x=query_length, xend=query_length, y=0, yend=real_time, color=kmer_color_label), size=kmer_size)

plot_e4 = plot_e4 + scale_x_continuous(name="k (query length)", breaks=seq(0, 31, 1)) + scale_y_continuous(name="runtime (ns)")
plot_e4 = plot_e4 + geom_hline(yintercept=0, color="grey")
plot_e4 = plot_e4 + ggtitle(label=title, subtitle=paste("text length =", text_length, "| mean of 150 benchmark cycles per k")) + theme(plot.title=element_text(face="bold"))

#1e5
########################################################################################
text_length = 1e5
data = data_raw[data_raw$text_length == text_length,]
data_kmer = data[data$name == "kmer_single_median",]
data_fm = data[data$name == "single_kmer_median",] #sic, typod the benchmark name

plot_e5 = ggplot()
plot_e5 = plot_e5 + geom_segment(data=data_fm, aes(x=query_length+fm_offset, xend=query_length+fm_offset, y=0, yend=real_time, color=fm_color_label), size=fm_size)
plot_e5 = plot_e5 + geom_segment(data=data_kmer, mapping=aes(x=query_length, xend=query_length, y=0, yend=real_time, color=kmer_color_label), size=kmer_size)

plot_e5 = plot_e5 + scale_x_continuous(name="k (query length)", breaks=seq(0, 31, 1)) + scale_y_continuous(name="runtime (ns)")
plot_e5 = plot_e5 + geom_hline(yintercept=0, color="grey")
plot_e5 = plot_e5 + ggtitle(label=title, subtitle=paste("text length =", text_length, "| mean of 150 benchmark cycles per k")) + theme(plot.title=element_text(face="bold"))

#1e6
########################################################################################
text_length = 1e6
data = data_raw[data_raw$text_length == text_length,]
data_kmer = data[data$name == "kmer_single_median",]
data_fm = data[data$name == "single_kmer_median",] #sic, typod the benchmark name

plot_e6 = ggplot()
plot_e6 = plot_e6 + geom_segment(data=data_fm, aes(x=query_length+fm_offset, xend=query_length+fm_offset, y=0, yend=real_time, color=fm_color_label), size=fm_size)
plot_e6 = plot_e6 + geom_segment(data=data_kmer, mapping=aes(x=query_length, xend=query_length, y=0, yend=real_time, color=kmer_color_label), size=kmer_size)

plot_e6 = plot_e6 + scale_x_continuous(name="k (query length)", breaks=seq(0, 31, 1)) + scale_y_continuous(name="runtime (ns)")
plot_e6 = plot_e6 + geom_hline(yintercept=0, color="grey")
plot_e6 = plot_e6 + ggtitle(label=title, subtitle=paste("text length =", text_length, "| mean of 150 benchmark cycles per k")) + theme(plot.title=element_text(face="bold"))

#1e7 #1e8 #1e9 #1e10
########################################################################################
plot_e6
ggsave("kmer_search_fm_vs_single.png", plot_e6)