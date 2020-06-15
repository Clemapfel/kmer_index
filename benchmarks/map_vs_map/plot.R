# Created by: clem
# Created on: 6/12/20

filter <- function(data, names, filter_percentage) {

    output <- data
    for (name in names) {

      temp <- data[data$name == name,]

      print(nrow(temp))

      mean <- mean(temp$real_time)
      stddev <- sd(temp$real_time)
      error <- qnorm(filter_percentage)*stddev/sqrt(nrow(temp))
      lower_bound <- mean - error
      upper_bound <- mean + error

      output <- output[
        (output$name == name & (output$real_time >= lower_bound & output$real_time <= upper_bound)) | output$name != name,]
    }

    return(output);
}

library(ggplot2)

setwd("/home/clem/Documents/Workspace/kmer_index/source/benchmarks/map_vs_map/")
data <- read.csv("2020-06-15_14-54-21.csv")

# filter out gbenchmark averages
data_at <- data[data$name == "da_at" | data$name == "std_at" | data$name == "robin_hood_at",]

plot_at <- qplot(name, real_time, data=data_at, fill=name, geom="boxplot")
plot_at = plot_at + ggtitle(label="unordered_map<size_t, std::vector<uint32_t>::at(size_t) performance",
                            subtitle="TODO")
plot_at = plot_at + theme(plot.title=element_text(face="bold"))
plot_at = plot_at + scale_fill_discrete(name = "name", labels = c("direct addressing", "robin_hood", "std"))

plot_at










data_insert <- data[data$name == "da_insert" | data$name == "std_insert" | data$name == "robin_hood_insert",]

plot_insert <- qplot(name, real_time, data=data_insert, fill=name, geom="boxplot")
plot_insert = plot_insert + ggtitle(label="unordered_map::insert(size_t, uint32_t) performance for different map implementations",
                                    subtitle="TODO")
plot_insert = plot_insert + theme(plot.title=element_text(face="bold"))
