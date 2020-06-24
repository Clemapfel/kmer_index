# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 6/19/20

library(ggplot2)

setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/")
data <- read.csv("final_1k_3k_2020-06-24T00-13-32+02-00.csv")
data_fm = data[data$name == "fm_median",]

#data <- read.csv("TODO:")
data_multi = data[data$name == "multi_kmer_median",]

plot <- ggplot()
#plot = plot + geom_area(mapping=aes(x=query_length, y=real_time, color="kmer_index"), data=data_multi, fill)
plot = plot + geom_segment(aes(x=query_length, xend=query_length, y=0, yend=data_multi$real_time, color="kmer_index"), data=data_multi, size=1)
plot = plot + geom_hline(yintercept=mean(data_multi$real_time))
#plot = plot + geom_hline(mapping=aes(yintercept=mean(data_fm$real_time, na.rm=TRUE), color="fm_index"))
#plot = plot + geom_hline(mapping=aes(yintercept=mean(data_multi$real_time, na.rm=TRUE), color="multi-kmer index"))

plot <- plot + geom_line(data = data_fm, mapping=aes(query_length, real_time, color="fm_index")) + geom_smooth(method="lm", data=data_fm, mapping=aes(query_length, real_time, color="fm_index"), size=2);
plot <- plot + coord_cartesian(ylim=c(min(data_multi$real_time), max(data_fm$real_time)), xlim=c(1000, 3000))
#plot <- plot + scale_y_continuous(breaks=seq(0, 6000, 500), name="runtime (ns)") + scale_x_continuous(breaks=seq(0, 200, 10), name="query length")
plot <- plot + ggtitle(label="exact search performance", subtitle="queries randomized each call, 10 benchmarks cycles per query length")
plot

plot_diff <- ggplot()
diff = data_multi$real_time - data_fm$real_time
plot_diff = plot_diff + geom_segment(mapping=aes(x=data_multi$query_length, y=0, xend=data_multi$query_length, yend=diff))
plot_diff = plot_diff + coord_cartesian(ylim=c(-2500, 2500))
plot_diff
