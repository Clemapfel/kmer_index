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
xlim = c(0, 45)

single_color = "blue"
multi_color = "navy"
fm_color = "red"

# MACRO
macro_line_size = 0.75

macro = ggplot()
macro = macro + geom_line(mapping=aes(x=query_length, y=real_time, color=single_color_label), data=data_single, size=macro_line_size)
#macro = macro + geom_line(mapping=aes(x=query_length, y=real_time, color=multi_color_label), data=data_multi, size=smooth_size)
macro = macro + scale_x_continuous(name="k (query length)", breaks=seq(1, 50, 1)) + scale_y_continuous(name="runtime (ns)", seq(0, 6e07, 1e07))
macro = macro + theme(legend.title=element_blank(), legend.position="none", plot.title=element_text(face="bold"))
macro = macro + ggtitle(label="search performance for kmer_index<10>", subtitle="queries randomized each call | 25 benchmark cycles per query length")

# inset plot zooming in on 15:20
data_zoomed <- data_single[data_single$query_length > 15 & data_single$query_length <= 20,]
x <- data_zoomed$query_length
y <- data_zoomed$real_time

inset_plot = ggplot()
inset_plot = inset_plot + geom_segment(mapping=aes(x=x, xend=x, y=0, yend=y, color=single_color_label), size=macro_line_size)

inset_plot <- inset_plot + theme(axis.title.y=element_blank(),
                                 axis.title.x=element_blank(),
                                 legend.position="none",
                                 plot.background=element_rect(color="black"))

inset_x = 25
inset_y = 6e+07
inset_width = 20
inset_height = 2.5e+07

inset_grob <- ggplotGrob(inset_plot)
macro = macro + annotation_custom(grob=inset_grob, xmin=inset_x, xmax=inset_x + inset_width, ymin=inset_y - inset_height, ymax=inset_y)
macro = macro + geom_segment(aes(x=16, y = 0, xend=inset_x, yend=inset_y, color="zoom", alpha=0.1))
macro = macro + geom_segment(aes(x=20, y = 0, xend=inset_x+inset_width, yend=inset_y - inset_height, colour="zoom", alpha=0.1))
macro = macro + geom_segment(aes(x=16, y = -0.1e+07, xend=16, yend=0.1e+07, color="zoom"), size=macro_line_size)
macro = macro + geom_segment(aes(x=19.8, y = -0.1e+07, xend=19.8, yend=0.1e+07, color="zoom"), size=macro_line_size)
macro = macro + scale_color_manual(values=c("red", "green", "blue", "black"), breaks=c(single_color_label, multi_color_label, fm_color_label, "zoom"))

macro


if (FALSE)
  {

multi_max = data_multi[data_multi$name=="multi_kmer_median" & data_multi$query_length==1,"real_time"]
micro = ggplot()
micro = micro + geom_line(mapping=aes(x=query_length, y=real_time, color=single_color_label), data=data_single, size=smooth_size-1)
micro = micro + geom_line(mapping=aes(x=query_length, y=real_time, color=multi_color_label), data=data_multi, size=smooth_size-1)
micro = micro + coord_cartesian(ylim=c(min(data_multi$real_time, multi_max), 40000), xlim=xlim)
micro = micro + scale_x_continuous(name="k (query length)", breaks=seq(1, 50, 2)) + scale_y_continuous(name="runtime (ns)", breaks=seq(0, 40000, 2000))
micro = micro + ggtitle(label="search performance for kmers", subtitle="single\t=\tkmer_index<10>\nmulti\t=\tkmer_index<5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29>")
micro = micro + theme(plot.title=element_text(face="bold"), legend.title=element_blank())


plot = ggplot()

#single
plot = plot + geom_line(mapping=aes(x=query_length, y=real_time, color=single_color_label), data=data_single, size=smooth_size-1)

#fm
#plot = plot + geom_line(mapping=aes(x=query_length, y=real_time, color=fm_color_label), data=data_fm)
#plot = plot + geom_smooth(method="lm", mapping=aes(x=query_length, y=real_time, color=fm_color_label), data=data_fm, size=smooth_size, fullrange=TRUE, se=FALSE)

#multi
plot = plot + geom_line(mapping=aes(x=query_length, y=real_time, color=multi_color_label), data=data_multi)
#plot = plot + geom_smooth(method="lm", mapping=aes(x=query_length, y=real_time, color=multi_color_label), data=data_multi[data_multi$query_length >= 4 & data_multi$query_length <= 34,], size=smooth_size, fullrange=TRUE, se=FALSE)

plot = plot + coord_cartesian(ylim=c(min(data_multi$real_time), 1e5), xlim=c(3, 47))
plot = plot + scale_x_continuous(name="k (query length)", breaks=seq(3, 50, 2)) + scale_y_continuous(name="runtime (ns)")
plot = plot + ggtitle(label="search performance for kmers", subtitle="fm\t\t=\tseqan3::fm-index\nsingle\t=\tkmer_index<10>\nmulti\t=\tkmer_index<5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29>")
plot = plot + theme(plot.title=element_text(face="bold"), legend.title=element_blank())
#plot = plot + scale_color_manual()
plot }