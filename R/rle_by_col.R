rle_by_col <- function(col) {
  run_encoded <- rle(as.integer(col))
  w <- which(run_encoded$values==1)
  if (w>1) {
    run_encoded$lengths[w+1] <- run_encoded$lengths[w+1]+run_encoded$lengths[w]-1
  } else {
    run_encoded$lengths[w+1] <- run_encoded$lengths[w]-1
    run_encoded$values[w+1] <- 0
  }
  run_encoded$values[w] <- run_encoded$lengths[w]
  run_encoded$lengths[w] <- 1
  return(inverse.rle(run_encoded))
}