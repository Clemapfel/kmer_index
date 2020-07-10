# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 6/19/20

library(ggplot2)
library(gridExtra)
library(grid)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/5_1000")
data = read.csv("1e6_experimental_2020-07-09T22-26-24+02-00.csv")
data = data[!grepl("stddev", data$name, fixed=TRUE) ,]

fm_color_label = "fm"
multi_color_label = "multi"

single_color = "deepskyblue2"
multi_color = single_color
fm_color = "lightcoral"

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
        output = c(output, 0)
    }
  }

  return(output);
}

get_diff = function(text_length) {

  data_fm = data[data$name == "fm_median" & data$text_length==text_length,]
  data_multi = data[data$name == "kmer_median" & data$text_length==text_length,]
  diff = data_fm$real_time - data_kmer$real_time
  #diff = diff * 100
  return(sign(diff) * log(sign(diff)*diff))
}

get_title = function(text_length) {

  return(ggtitle(label=paste("text length = ", text_length)))#, subtitle="fm\t\t=\tseqan3::fm-index\nmulti\t=\tkmer_index<5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>"))
}

stopifnot(FALSE);

ylim = c()
y_scale = scale_y_continuous(name = "log(speedup %)", breaks=seq(-10000000, 100000000, 1000))
x_scale = scale_x_continuous(name="query length", breaks=seq(0, 1000, 30))
theme = theme(plot.title=element_text(face="bold"), legend.position="none", axis.title.x=element_blank())
coord = coord_cartesian(ylim=c(-5,5))#min(diff), max(diff)))
color = scale_color_manual(name = "", values =c(fm_color, multi_color), labels = c("fm faster than kmer","kmer faster than fm"))

get_plot = function(text_length, col_a_name, col_b_name) {

  data_multi = data[data$name == col_a_name & data$text_length==text_length,]
  data_fm = data[data$name == col_b_name & data$text_length==text_length,]
  diff = speedup(data_multi$real_time, data_fm$real_time)
  #diff = sign(diff) * log(sign(diff) * diff * 100)
  query_length = seq(min(data_multi$query_length, na.rm=TRUE), max(data_multi$query_length), 1)

  speedup = diff
  plot = ggplot() + geom_segment(aes(x=query_length, xend=query_length, y=0, yend=speedup, color=ifelse(speedup>0, multi_color_label, fm_color_label)), size=2)
  plot = plot + ggtitle(label=paste("speedup(", col_a_name, col_b_name, ") | text length = ", text_length, sep="")) + y_scale + x_scale + theme + coord + color
}

print(get_plot(1e6, "kmer_median", "kmer_experimental_median"))

stopifnot(FALSE);

plot_1e4 = get_plot(1e4)
plot_1e6 = get_plot(1e6)
plot_1e8 = get_plot(1e8)

proxy = get_plot(1e4, 5) +  theme(legend.position="left", legend.text=element_text(size = 16))

plot = grid.arrange(plot_1e4, plot_1e6, plot_1e8, ncol=1,
                     top=textGrob(expression(bold("relative speedup (log-scaled): kmer vs. fm"))),
                     bottom=get_legend(proxy))#legendGrob(c("fm faster than kmer", "kmer faster than fm"), pch=c(15,15), gp=gpar(color=c(fm_color, multi_color_label))))

ggsave("runtime_diff_over_text_size.png", plot, width=30, height=20, units="cm")
plot


























if (FALSE) {

y_scale = scale_y_continuous(name = "speedup (%)", breaks=c(-1000, 1000, 10))
x_scale = scale_x_continuous(name="query length", breaks=seq(0, 1000, 30))
theme = theme(plot.title=element_text(face="bold"), legend.position="none", axis.title.x=element_blank())
coord = coord_cartesian(ylim=c(-100, 100))#min(diff), max(diff)))
color = scale_color_manual(name = "", values =c(fm_color, multi_color), labels = c("fm faster than kmer","kmer faster than fm"))


# 1e4
plot_1e4 = ggplot()
diff_1e4 = get_diff(1e4)
plot_1e4 = plot_1e4 + geom_segment(mapping=aes(x=query_length, xend=query_length, y=diff_1e4, yend=0, color=ifelse(diff_1e4<0, fm_color_label, multi_color_label)))
plot_1e4 = plot_1e4 + y_scale + x_scale + color + theme + coord + ggtitle(label=paste("text length = 1e+04"))
plot_1e4

# 1e6
plot_1e6 = ggplot()
diff_1e6 = get_diff(1e6)
plot_1e6 = plot_1e6 + geom_segment(mapping=aes(x=query_length, xend=query_length, y=diff_1e6, yend=0, color=ifelse(diff_1e6<0, fm_color_label, multi_color_label)))
plot_1e6 = plot_1e6 + y_scale + x_scale + color + theme + coord + get_title(1e6)
plot_1e6

# 1e8
plot_1e8 = ggplot()
diff_1e8 = get_diff(1e8)
plot_1e8 = plot_1e8 + geom_segment(mapping=aes(x=query_length, xend=query_length, y=diff_1e8, yend=0, color=ifelse(diff_1e8<0, fm_color_label, multi_color_label)))
plot_1e8 = plot_1e8 + y_scale + x_scale + color + theme + coord + get_title(1e8)
plot_1e8

proxy = ggplot() + geom_segment(mapping=aes(x=query_length, xend=query_length, y=diff_1e4, yend=0, color=ifelse(diff_1e4>0, fm_color_label, multi_color_label)))
proxy = proxy + color + theme(legend.title=element_blank(),  legend.text = element_text(size = 15))
#plot_1e4 = plot_1e4 + y_scale + x_scale + color + theme + coord + get_title(1e4)

# arrange
plot <- grid.arrange(plot_1e4, plot_1e6, plot_1e8, ncol=1,
                     top=textGrob(expression(bold("search runtime difference (kmer - fm) over query length"))),
                     bottom=get_legend(proxy))#legendGrob(c("fm faster than kmer", "kmer faster than fm"), pch=c(15,15), gp=gpar(color=c(fm_color, multi_color_label))))

plot
#ggsave("runtime_diff_over_text_size.png", plot, width=30, height=20, units="cm")
}

