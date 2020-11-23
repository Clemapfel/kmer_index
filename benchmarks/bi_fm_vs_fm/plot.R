# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 16.07.20

library(ggplot2)
library(ggpubr)
library(dplyr)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/bi_fm_vs_fm/")
data = read.csv("1e8_2020-07-16T15-44-04+02-00.csv")
#data = data[!grepl("stddev", data$name, fixed=TRUE),]

bi_fm = data[data$name=="bi_fm_median",]
fm = data[data$name=="fm_median",]

speedup = function(col_a, col_b) {

  output = c()
  for (i in seq_along(col_a)) {

    if (col_a[i] > col_b[i]) {
        output = c(output, 1 - col_a[i] / col_b[i])
    }
    else if (col_b[i] > col_a[i]) {
        output = c(output, -1*(1 - col_b[i] / col_a[i]))
    }
    else {
        output = c(output, 0);
    }
  }

  return(output);
}

size = 1

plot = ggplot() + geom_line(mapping=aes(x=fm$query_length, y=fm$real_time, color="fm-index"), size=3*size, alpha=0.7)
plot = plot + geom_line(aes(x=bi_fm$query_length, y=bi_fm$real_time, color="bi-fm-index"), size=1*size, alpha=1)
plot = plot + ggtitle(label="absolute runtime", subtitle="text size = 10^8 | 100 benchmark cycles per query length")
plot = plot + scale_y_continuous(name="runtime (ns)", breaks=seq(0, 1000000, 1000)) + scale_x_continuous(name="query length", breaks=seq(0,300,25))
plot = plot + theme(plot.title=element_text(face="bold"), legend.position="bottom", legend.title=element_blank())

diff = speedup(bi_fm$real_time, fm$real_time)*100
diff_plot = ggplot() + geom_segment(aes(x=fm$query_length, xend=fm$query_length, y=0, yend=diff, color=ifelse(diff>0,"bi_fm faster than fm", "fm faster than bi_fm")))
diff_plot = diff_plot + scale_y_continuous(name="speedup %", breaks=seq(-10,10,0.25)) + scale_x_continuous(name="query length", breaks=seq(0,300,25)) + coord_cartesian(ylim=c(-2,2))
diff_plot = diff_plot + ggtitle(label="speedup(bi_fm, fm)", subtitle=paste("mean speedup in [3,300] = ", round(mean(diff), digits=3), "%"))
diff_plot = diff_plot + theme(legend.position="bottom", legend.title=element_blank(), plot.title=element_text(face="bold"))
#diff_plot = diff_plot + scale_color_manual(labels=c("bi_fm faster than fm", "fm faster than bi_fm"))


plot = grid.arrange(plot, diff_plot, ncol=2, widths=c(10,20))
print(plot)
ggsave("bi_fm_vs_fm.png", plot, width=30, height=20, units="cm")



