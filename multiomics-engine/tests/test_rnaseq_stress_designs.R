#!/usr/bin/env Rscript
# Base-R regression of the prespecified RNA-seq stress generator.
source("multiomics-engine/rnaseq_stress_designs.R")
grid <- rnaseq_stress_grid()
stopifnot(nrow(grid)==4L,!anyDuplicated(grid$scenario),
          setequal(grid$batch_mode,c("balanced","unequal")),
          any(grid$abs_log2_fc<.5),
          any(grid$n_control<=6L),
          any(grid$n_control!=grid$n_treated),
          any(grid$n_alternatives==4L))
for(i in seq_len(nrow(grid))) {
  s <- grid[i,,drop=FALSE]
  x <- rnaseq_stress_generate(s,1L)
  y <- rnaseq_stress_generate(s,1L)
  z <- rnaseq_stress_generate(s,2L)
  stopifnot(identical(x$seed,y$seed),
    identical(x$x,y$x),identical(x$truth,y$truth),
    !identical(x$x,z$x),x$seed!=z$seed,
    nrow(x$x)==s$n_control+s$n_treated,ncol(x$x)==400L,
    identical(rownames(x$x),rownames(x$meta)),
    identical(colnames(x$x),names(x$truth)),
    !anyNA(x$x),all(is.finite(x$x)),all(x$x>=0),
    all(x$x==floor(x$x)),
    sum(x$truth!=0)==s$n_alternatives,
    sum(x$truth>0)==sum(x$truth<0),
    all(unique(abs(x$truth[x$truth!=0]))==s$abs_log2_fc),
    all(table(x$meta$condition,x$meta$batch)>=2L),
    qr(stats::model.matrix(~batch+condition,x$meta))$rank==3L,
    length(unique(x$meta$subject_id))==nrow(x$x),
    is.finite(x$observed_depth_ratio),
    is.finite(x$design_condition_number))
}
reject <- function(s,replicate=1L)
  inherits(try(rnaseq_stress_generate(s,replicate),silent=TRUE),"try-error")
bad <- grid[1,,drop=FALSE]
bad$abs_log2_fc <- 0
stopifnot(reject(bad))
bad <- grid[1,,drop=FALSE]
bad$n_alternatives <- 401L
stopifnot(reject(bad))
bad <- grid[1,,drop=FALSE]
bad$batch_mode <- "unequal"
stopifnot(reject(bad))
bad <- grid[3,,drop=FALSE]
bad$n_control <- 10L
stopifnot(reject(bad))
bad <- grid[1,,drop=FALSE]
bad$scenario <- "unregistered"
stopifnot(reject(bad))
stopifnot(reject(grid[1:2,,drop=FALSE]),
  reject(grid[1,,drop=FALSE],0L),
  reject(grid[1,,drop=FALSE],NA_integer_))
cat("RNA-seq stress grid: deterministic seeds, planted truth, independent samples and identifiable batch: PASS\n")
