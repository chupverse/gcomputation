.miboot_gc_binary <- function(formula, data, group, effect="ATE",
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
  
  p0 <- c()
  p1 <- c()
  OR <- c()
  delta <- c()
  ratio <- c()
  p0.unadj <- c()
  p1.unadj <- c()
  OR.unadj <- c()
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
    
    gc_res <- .gc_binary(formula = formula, data = current_imputed_data, group = group,
                          effect = effect, model = model, param.tune = param.tune,
                          cv = cv, boot.type = boot.type, boot.number = boot.number,
                          boot.tune = boot.tune, progress = FALSE, seed = gc_seed)
    
    datas[[i]] <- current_imputed_data
    calibrations[[i]] <- gc_res$calibration
    
    p0 <- c(p0, gc_res$p0)
    p1 <- c(p1, gc_res$p1)
    delta <- c(delta, gc_res$delta)
    ratio <- c(ratio, gc_res$ratio)
    OR <- c(OR, gc_res$OR)
    p0.unadj <- c(p0.unadj, gc_res$p0.unadj)
    p1.unadj <- c(p1.unadj, gc_res$p1.unadj)
    delta.unadj <- c(delta.unadj, gc_res$delta.unadj)
    ratio.unadj <- c(ratio.unadj, gc_res$ratio.unadj)
    OR.unadj <- c(OR.unadj, gc_res$OR.unadj)
    
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
  final_res$p0 <- p0
  final_res$p1 <- p1
  final_res$delta <- delta
  final_res$ratio <- ratio
  final_res$OR <- OR
  final_res$p0.unadj <- p0.unadj
  final_res$p1.unadj <- p1.unadj
  final_res$delta.unadj <- delta.unadj
  final_res$ratio.unadj <- ratio.unadj
  final_res$OR.unadj <- OR.unadj
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
