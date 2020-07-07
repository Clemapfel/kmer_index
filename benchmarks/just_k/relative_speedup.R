library(ggplot2)


# to 1e6
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/just_k/")
data_big = read.csv("new_result_complete.csv")#lncrna_complete.csv")
data_big = data_big[!grepl("stddev", data_big$name, fixed=TRUE) ,]
#data_big = data_big[!data_big$query_length<5, ]

speedup = function(col_a, col_b) {

  output = c()
  for (i in seq_along(col_a)) {

    if (col_a[i] > col_b[i]) {
        output = c(output, 1 - col_a[i] / col_b[i])
    }
    else if (col_b[i] > col_a[i]) {
        output = c(output, -1*(1 - col_b[i] / col_a[i]))
    }
  }

  return(output);
}

means = c()
ns = c()

for (text_length in sort(unique(data_big$text_length))) {

  data_cur = data_big[data_big$text_length == text_length,]
  fm = data_cur[data_cur$name=="fm_median",];
  kmer = data_cur[data_cur$name=="kmer_median",]
  print(paste("text length =", text_length));
  cur = speedup(kmer$real_time, fm$real_time)
  min_k = min(data_cur$query_length)
  for (s in cur) {
    if (s <= 0) {
      break
    }
    min_k = min_k + 1;
  }

  means = c(means, mean(cur)*100)
  ns = c(ns, text_length)
  print(paste("mean : ", round(mean(cur)*100, digits=2), "%", sep=""))
  print(paste("faster while k >", min_k))
}

plot = ggplot() + geom_line(aes(x=ns, y=means), alpha=0.5) + geom_smooth(aes(x=ns, y=means))



#####################

fm_color_label = "fm"
multi_color_label = "multi"

multi_color = "springgreen4"
fm_color = "red2"

get_speedup_plot = function(text_length) {

    fm = data_big[data_big$name == "fm_median" & data_big$text_length == text_length,]
    kmer = data_big[data_big$name == "kmer_median" & data_big$text_length == text_length,]
    percent = speedup(kmer$real_time, fm$real_time)*100

    plot = ggplot()
    plot = plot + geom_segment(aes(x=kmer$query_length, xend=kmer$query_length, y=0, yend=percent, color=ifelse(percent > 0, multi_color_label, fm_color_label)), size=4)
    plot = plot + scale_x_continuous(name="query length (k)", breaks=seq(1, 30, 1))
    plot = plot + scale_y_continuous(name="% speedup", breaks=seq(-1000, 1000, 5))
    plot = plot + ggtitle(label=paste("text size = ", text_length))
    plot = plot + theme(plot.title=element_text(face="bold"))
    plot = plot + geom_hline(yintercept=0, color="darkgrey")
    plot = plot + scale_color_manual(name = "", values=c(fm_color, multi_color), labels = c("fm faster than kmer","kmer faster than fm"))

    #plot = plot + coord_cartesian(ylim=c(-150,30))
    return(plot)
}

ggsave("relative_speedup.png", get_speedup_plot(1e8), width=30, height=15, units="cm")

get_diff_plot = function(text_length) {

    fm = data_big[data_big$name == "fm_median" & data_big$text_length == text_length,]
    kmer = data_big[data_big$name == "kmer_median" & data_big$text_length == text_length,]
    percent = kmer$real_time - fm$real_time

    plot = ggplot()
    plot = plot + geom_segment(aes(x=kmer$query_length, xend=kmer$query_length, y=0, yend=percent, color=ifelse(percent < 0, multi_color_label, fm_color_label)), size=4)
    plot = plot + scale_x_continuous(name="query length (k)", breaks=seq(1, 30, 1))
    plot = plot + scale_y_continuous(name="runtime difference (ns)", breaks=seq(-1000, 1000, 20))
    plot = plot + ggtitle(label="kmer - fm runtime", subtitle=paste("text size = ", text_length))
    plot = plot + theme(plot.title=element_text(face="bold"))
    plot = plot + geom_hline(yintercept=0, color="darkgrey")
    plot = plot + scale_color_manual(name = "", values=c(fm_color, multi_color), labels = c("fm faster than kmer","kmer faster than fm"))

    #plot = plot + coord_cartesian(ylim=c(-150,30))
    return(plot)
}

get_line_plot = function(text_length) {

    fm = data_big[data_big$name == "fm_median" & data_big$text_length == text_length,]
    kmer = data_big[data_big$name == "kmer_median" & data_big$text_length == text_length,]

    plot = ggplot()
    plot = plot + geom_line(aes(x=kmer$query_length, y=kmer$real_time, color=multi_color_label), size=2)
    plot = plot + geom_line(aes(x=fm$query_length, y=fm$real_time, color=fm_color_label), size=2)
    plot = plot + scale_x_continuous(name="query length (k)", breaks=seq(1, 30, 1))
    #plot = plot + scale_y_continuous(name="runtime (ns)", breaks=seq(-1000, 1000, 20))
    plot = plot + ggtitle(label=paste("text size = ", text_length))
    plot = plot + theme(plot.title=element_text(face="bold"))
    plot = plot + scale_color_manual(name = "", values=c(fm_color, multi_color), labels = c("fm index","kmer index"))

    #plot = plot + coord_cartesian(ylim=c(-150,30))
    return(plot)
}




