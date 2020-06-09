pick_best <- function(seq, modulos_vec) {
  output = c();
  for(i in seq) {
    output = c(output, min(i %% modulos_vec))
  }
  return(output)
}





