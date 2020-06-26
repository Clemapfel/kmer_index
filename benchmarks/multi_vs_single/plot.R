# Created by: clem
# Created on: 26.06.20

library(ggplot2)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_vs_single/")
data = read.csv("2020-06-26T18-07-30+02-00.csv")

data_fm = data[data$name=="fm_median",]
data_multi = data[data$name=="multi_kmer_median",]
data_single = data[data$name=="single_kmer_median",]

ks = c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29)
smooth_size = 2
fm_color_label = "fm"
single_color_label = "single"
multi_color_label = "multi"

plot = ggplot()

plot = plot + geom_line(mapping=aes(x=query_length, y=real_time, color=single_color_label), data=data_single, size=smooth_size-1)

plot = plot + geom_line(mapping=aes(x=query_length, y=real_time, color=fm_color_label), data=data_fm)
plot = plot + geom_smooth(method="lm", mapping=aes(x=query_length, y=real_time, color=fm_color_label), data=data_fm, size=smooth_size, fullrange=TRUE, se=FALSE)

plot = plot + geom_line(mapping=aes(x=query_length, y=real_time, color=multi_color_label), data=data_multi)
plot = plot + geom_smooth(method="lm", mapping=aes(x=query_length, y=real_time, color=multi_color_label), data=data_multi[data_multi$query_length >= 4 & data_multi$query_length <= 34,], size=smooth_size, fullrange=TRUE, se=FALSE)

plot = plot + coord_cartesian(ylim=c(min(data_multi$real_time), max(data_fm$real_time)+2000), xlim=c(3, 47))
plot = plot + scale_x_continuous(name="k (query length)", breaks=seq(3, 50, 2)) + scale_y_continuous(name="runtime (ns)")
plot = plot + ggtitle(label="search performance for kmers", subtitle="fm\t\t=\tseqan3::fm-index\nsingle\t=\tkmer_index<10>\nmulti\t=\tkmer_index<5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29>")
plot = plot + theme(plot.title=element_text(face="bold"), legend.title=element_blank())
#plot = plot + scale_color_manual()
plot