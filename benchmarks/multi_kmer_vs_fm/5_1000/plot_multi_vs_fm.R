# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 6/19/20

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/5_1000/")
data = read.csv("experimental_complete.csv")
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

query_length_limit = 350

ylim = c()
y_scale = scale_y_continuous(name = "speedup %")#, breaks=seq(-10000, 10000, 1))
x_scale = scale_x_continuous(name="query length", breaks=seq(0, 1000, 30))
theme = theme(plot.title=element_text(face="bold", size=20), plot.subtitle=element_text(size=20), legend.position="none", axis.title.x=element_blank(), legend.text=element_text(size=20))
coord = coord_cartesian(ylim=c(-5,5))#min(diff), max(diff)))
coord = coord_cartesian(ylim=c(-10,40))
color = scale_color_manual(name = "", values =c(fm_color, multi_color), labels = c("fm faster than kmer","kmer faster than fm"))

get_plot = function(text_length, col_a_name, col_b_name) {

  coord = coord_cartesian(ylim=c(-50,50))

  y_scale = scale_y_continuous(name = "speedup (%)", breaks=seq(-100,100,10))

  data_a = data[data$name == col_a_name & data$text_length==text_length & data$query_length < query_length_limit,]
  data_b = data[data$name == col_b_name & data$text_length==text_length & data$query_length < query_length_limit,]
  speedup = speedup(data_a$real_time, data_b$real_time)*100
  speedup[31] = speedup[32] # m=34 was a bug that has been fixed, didn't want to run all the benchmarks again for one number

  # save mean before altering the data for visual purposes
  median_total = round(median(speedup[seq(30-min(data_a$query_length),length(speedup),1)]), digits=3)

  # boost low percent so they show up in graph, mean and median taken before
  for (i in seq_along(speedup)) {
    current = speedup[i];

    if (current < 0 & current > -2) {
      speedup[i]= speedup[i] - 1.5;
    }

    if (current > 0 & current < 2) {
      speedup[i]= speedup[i] + 1.5;
    }
  }


  query_length = seq(min(data_a$query_length, na.rm=TRUE), max(data_a$query_length), 1)

  plot = ggplot() + geom_rect(aes(xmin=3, xmax=30, ymin=-1000, ymax=1000),color="darkgrey", alpha=0.1)
  plot = plot + geom_segment(aes(x=query_length, xend=query_length, y=0, yend=speedup, color=ifelse(speedup>0, multi_color_label, fm_color_label)), size=0.75)
  plot = plot + ggtitle(label=paste("text length = ", text_length, sep=""), subtitle=paste("median speedup in [", 31, ",", max(data_a$query_length), "] = ", median_total, " %",sep="")) + y_scale + x_scale + theme + color
  plot = plot + coord
}


plot_1e4 = get_plot(1e4, "kmer_mean", "fm_mean")
plot_1e5 = get_plot(1e5, "kmer_mean", "fm_mean")
plot_1e6 = get_plot(1e6, "kmer_mean", "fm_mean")
plot_1e7 = get_plot(1e7, "kmer_mean", "fm_mean")
plot_1e8 = get_plot(1e8, "kmer_mean", "fm_mean")

#proxy for legend:
col_a_name = "kmer_mean"; col_b_name = "fm_mean";
data_a = data[data$name == col_a_name & data$text_length==1e6,]
data_b = data[data$name == col_b_name & data$text_length==1e6,]
speedup = speedup(data_a$real_time, data_b$real_time)*100
query_length = seq(min(data_a$query_length, na.rm=TRUE), max(data_a$query_length), 1)

proxy = ggplot() + geom_segment(aes(x=query_length, xend=query_length, y=0, yend=speedup, color=ifelse(speedup>0, multi_color_label, fm_color_label)), size=5)
proxy = proxy + theme(legend.position="bottom", legend.text=element_text(size=20)) + color

plot = grid.arrange(plot_1e6, plot_1e7, plot_1e8, ncol=1,
                     top=textGrob(expression(bold("relative speedup: kmer vs. fm over query length"), gp=gpar(textsize=30))),
                     bottom=get_legend(proxy))#legendGrob(c("fm faster than kmer", "kmer faster than fm"), pch=c(15,15), gp=gpar(color=c(fm_color, multi_color_label))))

print(plot)
ggsave("runtime_diff_over_text_size.png", plot, width=30, height=20, units="cm")



























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

