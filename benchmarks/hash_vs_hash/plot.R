# Created by: clem
# Created on: 6/8/20

# plot runtime for pows
library(ggplot2)
library(gridExtra)
library(grid)

filter_extremes <- function(data, interval=0.975) {

    # hash for
    data_for <- data[data$name == "hash_for",]
    mean_for <- mean(data_for$real_time)
    stddev_for <- sd(data_for$real_time)
    error_for <- qnorm(interval)*stddev_for/sqrt(nrow(data_for))
    lower_for <- mean_for - error_for
    upper_for <- mean_for + error_for

    # hash fold
    data_fold <- data[data$name == "hash_fold",]
    mean_fold <- mean(data_fold$real_time)
    stddev_fold <- sd(data_fold$real_time)
    error_fold <- qnorm(interval)*stddev_fold/sqrt(nrow(data_fold))
    lower_fold <- mean_fold - error_fold
    upper_fold <- mean_fold + error_fold

    return(data[(data$name == "hash_for" & data$real_time < upper_for & data$real_time > lower_for)
               |(data$name == "hash_fold" & data$real_time < upper_fold & data$real_time > lower_fold),])
}

round_to_nice <- function(x) {

    remainder = x %% 5
    return(x+10-remainder-10)
}

create_plot <- function(data, k, yscale=30, use_lab = FALSE) {

    temp <- data[data$k == k,]
    temp <- filter_extremes(temp,interval)
    min = round_to_nice(min(temp$real_time))
    y_seq = seq(min, min+yscale, 5)

    plot <- qplot(name, real_time, data=temp, fill=name, geom="boxplot", main=paste("k = ", k), xlab="")
    plot <- plot + theme(legend.position="none") + stat_summary(fun="mean", geom="point", shape=mean_point_shape, size=mean_point_size)

    if (use_lab) {
        plot = plot + scale_y_continuous(name = "real time (ns)", labels=y_seq, breaks=y_seq, limit=c(min, min+yscale))
        first = FALSE
    }
    else {
        plot = plot + scale_y_continuous(name="", labels=y_seq, breaks=y_seq, limit=c(min, min+yscale))
    }

    return(plot);

}

setwd("/home/clem/Workspace/kmer_index/source/benchmarks/hash_vs_hash/")
data_raw <- read.csv("2020-06-11_16-15-36.csv", stringsAsFactors = FALSE)

# filter out gbenchmark averages
data <- data_raw[data_raw$name=="hash_fold" | data_raw$name=="hash_for", ]

interval = 0.999
mean_point_size = 6
mean_point_shape = 4

yscale = 35;

plot5 <- create_plot(data, 5, yscale, TRUE)
plot10 <- create_plot(data, 10, yscale)
plot15 <- create_plot(data, 15, yscale)
plot20 <- create_plot(data, 20, yscale)
plot25 <- create_plot(data, 25, yscale)

plot <- grid.arrange(plot5, plot10, plot15, plot20, plot25, ncol=5,
                     top=textGrob(expression("hash performance for queries of length k\n(100 benchmark cycles per hash implementation and k, queries randomized each call)")),
                     right=legendGrob(c("mean"), pch=c(mean_point_shape)))

ggsave("hash_vs_hash.png", plot, width=30, height=15, units="cm")





