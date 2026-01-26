.miboot_gc_count <- function(formula, data, group, effect="ATE",
                                  model, param.tune=NULL, cv=10,
                                  boot.type="bcv", boot.number=500, boot.tune=FALSE,
                                  progress=TRUE, seed=NULL, m=5, ...) {
  
  cl <- match.call()
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (any(is.na(data))) {
    nmiss_initial <- nrow(data) - nrow(na.omit(data))
  } else {
    nmiss_initial <- 0
  }
  
  mice_args <- list(...)
  mice_args$data <- data
  mice_args$m <- m
  mice_args$printFlag <- FALSE
  
  imp_object <- do.call(mice::mice, mice_args)
  
  c0 <- c()
  c1 <- c()
  
  delta <- c()
  ratio <- c()
  c0.unadj <- c()
  c1.unadj <- c()
  
  delta.unadj <- c()
  ratio.unadj <- c()
  lambdas <- c()
  alphas <- c()
  datas <- list()
  calibrations <- list()
  
  final_res <- NULL
  
  if (progress) {
    pb <- utils::txtProgressBar(min = 0, max = m, style = 3, width = 50, char = "=")
    ip <- 0
  }
  
  for (i in 1:m) {
    if (progress) {
      ip <- ip + 1
      utils::setTxtProgressBar(pb, ip)
    }
    
    current_imputed_data <- mice::complete(imp_object, i)
    
    gc_seed <- if (!is.null(seed)) (seed + i) else NULL
    
    gc_res <- .gc_count(formula = formula, data = current_imputed_data, group = group,
                             effect = effect, model = model, param.tune = param.tune,
                             cv = cv, boot.type = boot.type, boot.number = boot.number,
                             boot.tune = boot.tune, progress = FALSE, seed = gc_seed)
    
    datas[[i]] <- current_imputed_data
    calibrations[[i]] <- gc_res$calibration
    
    c0 <- c(c0, gc_res$adjusted.results$c0)
    c1 <- c(c1, gc_res$adjusted.results$c1)
    delta <- c(delta, gc_res$adjusted.results$delta)
    ratio <- c(ratio, gc_res$adjusted.results$ratio)
    c0.unadj <- c(c0.unadj, gc_res$unadjusted.results$c0)
    c1.unadj <- c(c1.unadj, gc_res$unadjusted.results$c1)
    delta.unadj <- c(delta.unadj, gc_res$unadjusted.results$delta)
    ratio.unadj <- c(ratio.unadj, gc_res$unadjusted.results$ratio)
    
    if (model %in% c("lasso", "ridge")) {
      if (!is.null(gc_res$tuning.parameters$lambda)) {
        lambdas <- c(lambdas, gc_res$tuning.parameters$lambda)
      }
    } else if (model == "elasticnet") {
      if (!is.null(gc_res$tuning.parameters$lambda)) {
        lambdas <- c(lambdas, gc_res$tuning.parameters$lambda)
      }
      if (!is.null(gc_res$tuning.parameters$alpha)) {
        alphas <- c(alphas, gc_res$tuning.parameters$alpha)
      }
    }
    
    if (is.null(final_res)) {
      final_res <- gc_res
    }
  }
  
  if (progress) {
    close(pb)
  }
  
  
  final_res$data <- datas
  final_res$calibration <- calibrations
  
  final_res$adjusted.results <- data.frame(c1 = c1, c0 = c0, delta = delta, ratio = ratio)
  final_res$unadjusted.results <- data.frame(c1 = c1.unadj, c0 = c0.unadj, delta = delta.unadj, ratio = ratio.unadj)
  
  final_res$boot.number <- m * boot.number
  final_res$seed <- seed
  final_res$call <- cl
  final_res$m <- m
  final_res$initial.data <- data
  final_res$nimput <- nmiss_initial
  
  if (model %in% c("lasso", "ridge")) {
    final_res$tuning.parameters=list(lambda=lambdas)
  } else if (model == "elasticnet") {
    final_res$tuning.parameters=list(alpha=alphas, lambda=lambdas)
  }
  
  
  return(final_res)
}
