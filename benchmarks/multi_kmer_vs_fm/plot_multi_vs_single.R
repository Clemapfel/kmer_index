# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 6/17/20

library(ggplot2)

setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/5_1000/")
data_fm_single <- read.csv("experimental_complete.csv")

data_single = data_fm_single[data_fm_single$name == "single_kmer_median",]
data_fm = data_fm_single[data_fm_single$name == "fm_median",]

q_lines = ggplot() + geom_line(data=data_single, mapping=aes(x=query_length, y=real_time, color="kmer_index<10>"))
q_lines = q_lines + geom_line(data=data_multi, mapping=aes(x=query_length, y=real_time, color="kmer_index<10, 13, 15, 17>"))
#q_lines = q_lines + geom_smooth(data=data_fm, mapping=aes(x=query_length, y=real_time, color="seqan3::fm-index"))
y_lim = max(data_multi$real_time)
q_lines = q_lines + scale_x_continuous(breaks=seq(0, 200, 10)) + coord_cartesian(ylim=c(0, y_lim), xlim=c(5, 150))
q_lines
######################
if (FALSE)
  {

# if a is faster -> positive value
to_ratio <- function(a, b)
{
  stopifnot(length(a) == length(b))

  output <- c();
  for (i in 1:length(a)) {

    if (identical(a[i], NaN) | identical(b[i], NaN) | is.na(a[i]) | is.na(b[i])) {
      output = c(output, NaN)
      next;
    }

    print(paste(a[i], " ", b[i]))

    if (a[i] > b[i]) {
      output = c(output, (-1*a[i]/b[i])-1)
    }
    else if (a[i] < b[i]) {
      output = c(output, (1*b[i]/a[i])+1)
    }
    else
      output = c(output, 1)
  }

  return(output);
}

multi = data_multi$real_time
single = data_single$real_time

y <- to_ratio(multi,single)#data_multi$real_time / data_single$real_time;
x <- data_single$query_length

q_ratio <- qplot(x=x, xend=x, y=0, yend=y, geom="segment", xlab="query length", main="multi-kmer_index<10, 13, 15, 17> speedup compared to single_kmer_index<10>")
q_ratio <- q_ratio + scale_y_continuous(name="y times speedup", breaks=seq(-100, 100, 0.5));
q_ratio <- q_ratio + scale_x_continuous(name="query length", breaks=seq(0, 200, 5))
q_ratio <- q_ratio + geom_hline(yintercept=0, color="grey") + coord_cartesian(ylim=c(-100, 100), xlim=c(5, 200))

}
