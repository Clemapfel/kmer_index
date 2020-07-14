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
    else {
        output = c(output, 0);
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
  print(paste("max : ", round(max(cur)*100, digits=2), " at " , match(max(cur), cur)+min(data_cur$query_length)-1))
  print(paste("min : ", round(min(cur)*100, digits=2), " at " , match(min(cur), cur)+min(data_cur$query_length)-1))
  print(paste("faster while k >", min_k))
}

plot = ggplot() + geom_line(aes(x=ns, y=means), alpha=0.5) + geom_smooth(aes(x=ns, y=means))



#####################

fm_color_label = "fm"
multi_color_label = "multi"

single_color = "deepskyblue2"
multi_color = single_color
fm_color = "lightcoral"

get_line_plot = function(text_length) {

    fm = data_big[data_big$name == "fm_median" & data_big$text_length == text_length,]
    kmer = data_big[data_big$name == "kmer_median" & data_big$text_length == text_length,]

    linesize = 1.5
  alpha = 0.6

    plot = ggplot()
    plot = plot + geom_line(aes(x=fm$query_length, y=fm$real_time, color=fm_color_label), size=linesize, alpha=alpha)
    plot = plot + geom_line(aes(x=kmer$query_length, y=kmer$real_time, color=multi_color_label), size=linesize, alpha=alpha)
    plot = plot + scale_x_continuous(name="query length (k)", breaks=seq(1, 30, 2))
    plot = plot + scale_y_continuous(name="runtime (ns)", breaks=seq(-100000, +100000, 100))
    plot = plot + ggtitle(label="absolute runtime")
    plot = plot + theme(plot.title=element_text(face="bold", size=10), legend.position="bottom", panel.background =element_rect(fill="gray90"))
    plot = plot + scale_color_manual(name = "", values=c("coral", "skyblue2"), labels = c("fm index","kmer index"))

    #plot = plot + coord_cartesian(ylim=c(-150,30))
    return(plot)
}

get_diff_plot = function(text_length) {

    fm = data_big[data_big$name == "fm_median" & data_big$text_length == text_length,]
    kmer = data_big[data_big$name == "kmer_median" & data_big$text_length == text_length,]

    linesize = 1.5
    alpha = 0.8
    diff = kmer$real_time - fm$real_time

    plot = ggplot()
    plot = plot + geom_segment(aes(x=fm$query_length, xend=fm$query_length, y=0, yend=diff, color=ifelse(diff > 0, multi_color_label, fm_color_label)), size=linesize, alpha=alpha)
    plot = plot + scale_x_continuous(name="query length (k)", breaks=seq(1, 30, 2))
    plot = plot + scale_y_continuous(name="runtime (ns)", breaks=seq(-100000, +100000, 100))
    plot = plot + ggtitle(label="(kmer - fm) absolute runtime difference")
    plot = plot + theme(legend.position="bottom", panel.background =element_rect(fill="gray90"))
    plot = plot + scale_color_manual(name = "", values=c("coral", "skyblue2"), labels = c("fm index","kmer index"))

    #plot = plot + coord_cartesian(ylim=c(-150,30))
    return(plot)
}

get_speedup_plot = function(text_length) {

    fm = data_big[data_big$name == "fm_median" & data_big$text_length == text_length,]
    kmer = data_big[data_big$name == "kmer_median" & data_big$text_length == text_length,]
    percent = speedup(kmer$real_time, fm$real_time)*100

    plot = ggplot()
    plot = plot + geom_hline(yintercept=0, color="darkgrey", size=1)
    plot = plot + geom_segment(aes(x=kmer$query_length, xend=kmer$query_length, y=0, yend=percent, color=ifelse(percent < 0, fm_color_label, multi_color_label)), size=4)
    plot = plot + scale_x_continuous(name="query length (k)", breaks=seq(1, 30, 1))
    plot = plot + scale_y_continuous(name="% speedup", breaks=seq(-1000, 1000, 5))
    plot = plot + ggtitle(label="speedup(kmer, fm)", subtitle=paste("text size = ", text_length))
    plot = plot + theme(plot.title=element_text(face="bold"), legend.position="bottom")
    plot = plot + scale_color_manual(name = "", values=c(fm_color, multi_color), labels = c("fm faster than kmer","kmer faster than fm"))


    inset_grob = ggplotGrob(get_line_plot(text_length))
    plot = plot + annotation_custom(grob=inset_grob, xmin=17, xmax=30, ymin=max(percent)-40, ymax=max(percent))

    #plot = plot + coord_cartesian(ylim=c(-150,30))
    return(plot)
}

print(get_speedup_plot(1e8))
ggsave("relative_speedup.png", get_speedup_plot(1e8), width=30, height=20, units="cm")






