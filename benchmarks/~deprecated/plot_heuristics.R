data = read.csv(paste(getwd(), "/source/benchmarks/csv_dumps/heuristics.csv", sep=""))
heuristic = data[data$name== "heuristic",]
no_heuristic = data[data$name=="no_heuristic",]
plot(heuristic$k, heuristic$iterations, type="l", col="blue")
lines(no_heuristic$k, no_heuristic$iterations, type="l", col="red")
legend(x=8, y=2500, legend=c("with heuristic", "without heuristic"), col=c("blue", "red"), lty=1:2)
