# Title     : TODO
# Objective : TODO
# Created by: clem
# Created on: 6/19/20

library(ggplot2)
setwd("/home/clem/Workspace/kmer_index/source/benchmarks/multi_kmer_vs_fm/")
data <- read.csv("new_2020-06-29T19-30-50+02-00.csv")


fm_color_label = "fm"
single_color_label = "single"
multi_color_label = "multi"
xlim = c(0, 45)

multi_color = "springgreen4"
fm_color = "red2"

all_primes= c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,
199,211,223,227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311,
313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421,431,
433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,
563,569,571,577,587,593,599,601,607,613,617,619,631,641,643,647,653,659,661,
673,677,683,691,701,709,719,727,733,739,743,751,757,761,769,773,787,797,809,
811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,
941,947,953,967,971,977,983,991,997)


for (text_length in c(1e6)){#}, 1e9)) {
data_fm = data[data$name == "fm_median" & data$text_length==text_length,]
data_multi = data[data$name == "multi_kmer_median" & data$text_length==text_length,]
plot_diff <- ggplot()
diff = data_multi$real_time - data_fm$real_time
diff = sign(diff) * log(sign(diff)*diff)

all_ks = c(5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31)

covered = c();

for (i in 1:1000) {
  is_covered = FALSE;

  for (k in all_ks) {
    if (i%%k == 0) {
      is_covered=TRUE
    }
  }

  if (!is_covered){
      covered = c(covered, i)
  }
}

  bad_runtime = c()
  for (i in 1:length(diff)) {
    if (diff[i] > 0) {
      bad_runtime=c(bad_runtime, i)
    }
  }

  #plot_diff = plot_diff + geom_segment(aes(x=covered, xend=covered, y=-50, yend=+50), color="darkgrey")



#plot_diff = plot_diff + geom_segment(mapping=aes(x=all_primes, xend=all_primes, y=-100, yend=+100, color="primes"), size=1, alpha=0.3)
plot_diff = plot_diff + geom_segment(mapping=aes(x=data_multi$query_length, xend=data_multi$query_length, y=diff, yend=0, color=ifelse(diff>0, fm_color_label, multi_color_label)))

#plot_diff = plot_diff + geom_segment(mapping=aes(x=data_multi$query_length, y=0, xend=data_multi$query_length, yend=diff))
plot_diff = plot_diff + coord_cartesian(ylim=c(min(diff), max(diff)))
plot_diff = plot_diff + scale_y_continuous(name = "sign(delta) * log(sign(delta))") + scale_x_continuous(name="query length", breaks=bad_runtime)
plot_diff = plot_diff + ggtitle(label=paste("search time: (kmer - fm) | text length = ", text_length), subtitle="fm\t\t=\tseqan3::fm-index\nmulti\t=\tkmer_index<5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31>") + theme(plot.title=element_text(face="bold"))
plot_diff = plot_diff + scale_color_manual(name = "", values =c(fm_color, multi_color), labels = c("fm faster than kmer","kmer faster than fm"))
plot_diff
ggsave(paste("kmer_-_fm_", text_length, ".png"), plot_diff, width=100, height=15, units="cm")

}
plot_diff



