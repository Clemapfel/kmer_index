# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 10.07.20

library(ggplot2)
library(gridExtra)
library(grid)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/5_1000")
data = read.csv("1e6_experimental_2020-07-10T18-47-34+02-00.csv")
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

y_scale = scale_y_continuous(breaks=seq(-100, 100, 1))
x_scale = scale_x_continuous(breaks=seq(0, 1000, 30))

data_fm = data[data$name == "fm_median",]
data_reg = data[data$name == "kmer_median",]
data_exp = data[data$name == "kmer_experimental_median",]

plot_reg_fm = ggplot() + geom_segment(aes(x=data_fm$query_length, xend=data_fm$query_length, y=0, yend=speedup(data_reg$real_time, data_fm$real_time)*100, color="fm_reg"))
plot_reg_fm = plot_reg_fm + ggtitle(label="speedup(kmer_reg, fm)") + y_scale + x_scale
plot_reg_fm = plot_reg_fm + coord_cartesian(ylim=c(-10, 10))


plot_exp_fm = ggplot() + geom_segment(aes(x=data_fm$query_length, xend=data_fm$query_length, y=0, yend=speedup(data_exp$real_time, data_fm$real_time)*100, color="fm_exp"))
plot_exp_fm = plot_exp_fm + ggtitle(label="speedup(kmer_exp, fm)") + y_scale + x_scale
plot_exp_fm = plot_exp_fm + coord_cartesian(ylim=c(-10, 10))

plot_exp_reg = ggplot() + geom_segment(aes(x=data_fm$query_length, xend=data_fm$query_length, y=0, yend=speedup(data_exp$real_time, data_reg$real_time)*100, color="reg_exp"))
plot_exp_reg = plot_exp_reg + ggtitle(label="speedup(kmer_exp, kmer_reg)") + x_scale
plot_exp_reg = plot_exp_reg + coord_cartesian(ylim=c(-10, 300))
plot_exp_reg

