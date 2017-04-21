args <- commandArgs(TRUE)

get_num <- function(fpath){
  eq <- c()
  lf1 <- list.files(fpath)
  lf2 <- grep("ge", lf1)
  for(i in 1:length(lf2)){
    ge1 <- readRDS(paste(fpath, lf1[lf2][[nums[i]]], sep = ""))
    if(any(is.na(ge1))){eq[i] <- FALSE;next}
    if(any(is.na(ge1$eqst))){eq[i] <- FALSE;next}
    
    eq[i] <- TRUE
  }
  return(eq)
}


sum(get_num(args))