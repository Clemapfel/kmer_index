# Created by: clem
# Created on: 6/8/20

library(ggplot2)

setwd("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/kmer_vs_fm_over_query_size/")
data <- read.csv("k10_2020-06-04_17-14-01.csv")
data_median <- data[data$name=="kmer_search_median",]

plot <- qplot(data_median$query_length,
              data_median$real_time,
              geom = "line",
              xlab = "query length",
              ylab = "runtime (ns)")

plot = plot + ggtitle(label="kmer index exact search runtime per query length for k = 10, text length = 1000000",
                      subtitle = "(5 benchmark cycles per query length, search queries randomized each call)")
plot = plot + theme(plot.title=element_text(face="bold"))

# add horizontal lines for each k
for (i in 1:6) {

  plot = plot + geom_vline(xintercept=10*i, color="darkgrey", linetype="dotted", size=1)
}

plot <- plot + scale_x_continuous(breaks=append(9:25, seq(30, 60, 10)))

# inset plot zooming in on 15:20
data_zoomed <- data_median[data_median$query_length > 15 & data_median$query_length <= 20,]
x <- data_zoomed$query_length
y <- data_zoomed$real_time

inset_plot <- qplot(x=x, xend=x, y=0, yend=y, data=data_zoomed,
                   geom="segment", xlab="query length")
inset_plot <- inset_plot + theme(axis.title.y=element_blank(),
                                plot.background=element_rect(color="black"))
# inset top left
inset_x = 32
inset_y = 5e+07
inset_width = 25
inset_height = 2.5e+07

# add lines for zoom
inset_grob <- ggplotGrob(inset_plot)
plot = plot + annotation_custom(grob=inset_grob, xmin=inset_x, xmax=inset_x + inset_width, ymin=inset_y - inset_height, ymax=inset_y)
plot = plot + geom_segment(aes(x=16, y = 0, xend=inset_x, yend=inset_y, color="black", alpha=0.1))
plot = plot + geom_segment(aes(x=20, y = 0, xend=inset_x+inset_width, yend=inset_y - inset_height, colour="black", alpha=0.1))
plot = plot + geom_segment(aes(x=16, y = -0.1e+07, xend=16, yend=0.1e+07, color="black"))
plot = plot + geom_segment(aes(x=19.8, y = -0.1e+07, xend=19.8, yend=0.1e+07, color="black"))
plot = plot + theme(legend.position="none")

ggsave("k10_search_over_query.png", plot, width=30, height=15, units="cm")




