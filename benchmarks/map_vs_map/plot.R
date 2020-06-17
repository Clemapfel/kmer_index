# Created by: clem
# Created on: 6/12/20

library(ggplot2)
library(tidyr)

setwd("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/map_vs_map/")
data <- read.csv("2020-06-15_17-54-44.csv")


plot = ggplot(data=data)
plot = plot + geom_line(mapping=aes(size, real_time, color=name), data=data[data$name=="std_median",], size=2)
plot = plot + geom_line(mapping=aes(size, real_time, color=name), data=data[data$name=="boost_median",], size=2)
plot = plot + geom_line(mapping=aes(size, real_time, color=name), data=data[data$name=="abseil_median",], size=2)
plot = plot + geom_line(mapping=aes(size, real_time, color=name), data=data[data$name=="robin_hood_median",], size=2)
plot = plot + ggtitle(label="unordered_map<size_t, std::vector<uint32_t>>::at performance for different implementations", subtitle = "runtime over size of the map, sampled at n*50000\nqueries randomized each call, 20 benchmark cycles per map size")
plot = plot + xlab("# entries in map") + ylab("runtime (ns)") + xlim(c(0, 1.5*1e+06))
plot = plot + theme(plot.title=element_text(face="bold"))
labels = c("absl::", "boost::", "robin_hood::", "std::")
colors = c("darkorange1", "lightskyblue", "red", "blue")
plot = plot + scale_color_manual(name="namespace", labels=labels, values=colors)
plot
ggsave("map_vs_map.png", plot, width=30, height=15, units="cm")
