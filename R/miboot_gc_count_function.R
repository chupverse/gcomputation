.miboot_gc_count <- function(formula, data, group, effect="ATE",
                                  model, param.tune=NULL, cv=10,
                                  boot.type="bcv", boot.number=500, boot.tune=FALSE,
                                  progress=TRUE, seed=NULL, m=5, ...) {
  
  cl <- match.call()
  
  if (!is.null(seed)) {
    set.seed(seed)
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
    
    c0 <- c(c0, gc_res$c0)
    c1 <- c(c1, gc_res$c1)
    delta <- c(delta, gc_res$delta)
    ratio <- c(ratio, gc_res$ratio)
    c0.unadj <- c(c0.unadj, gc_res$c0.unadj)
    c1.unadj <- c(c1.unadj, gc_res$c1.unadj)
    delta.unadj <- c(delta.unadj, gc_res$delta.unadj)
    ratio.unadj <- c(ratio.unadj, gc_res$ratio.unadj)
    
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
  final_res$c0 <- c0
  final_res$c1 <- c1
  final_res$delta <- delta
  final_res$ratio <- ratio
  final_res$c0.unadj <- c0.unadj
  final_res$c1.unadj <- c1.unadj
  final_res$delta.unadj <- delta.unadj
  final_res$ratio.unadj <- ratio.unadj
  final_res$boot.number <- m * boot.number
  final_res$seed <- seed
  final_res$call <- cl
  final_res$m <- m
  final_res$initial.data <- data
  
  if (model %in% c("lasso", "ridge")) {
    final_res$tuning.parameters=list(lambda=lambdas)
  } else if (model == "elasticnet") {
    final_res$tuning.parameters=list(alpha=alphas, lambda=lambdas)
  }
  
  
  return(final_res)
}
