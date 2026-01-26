.miboot_gc_binary <- function(formula, data, group, effect="ATE",
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
  formulas <- list()
  
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
    
    p0 <- c(p0, gc_res$adjusted.results$p0)
    p1 <- c(p1, gc_res$adjusted.results$p1)
    delta <- c(delta, gc_res$adjusted.results$delta)
    ratio <- c(ratio, gc_res$adjusted.results$ratio)
    OR <- c(OR, gc_res$adjusted.results$OR)
    p0.unadj <- c(p0.unadj, gc_res$unadjusted.results$p0)
    p1.unadj <- c(p1.unadj, gc_res$unadjusted.results$p1)
    delta.unadj <- c(delta.unadj, gc_res$unadjusted.results$delta)
    ratio.unadj <- c(ratio.unadj, gc_res$unadjusted.results$ratio)
    OR.unadj <- c(OR.unadj, gc_res$unadjusted.results$OR)
    
    if (model %in% c("all","aic","bic")) {
      formulas[[i]] <- gc_res$tuning.parameters
    }
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
  final_res$adjusted.results <- list(p1 = p1, p0 = p0, delta = delta, ratio = ratio, OR = OR)
  final_res$unadjusted.results <- list(p1 = p1.unadj, p0 = p0.unadj, delta = delta.unadj, ratio = ratio.unadj, OR = OR.unadj)
  final_res$boot.number <- m * boot.number
  final_res$seed <- seed
  final_res$call <- cl
  final_res$m <- m
  final_res$initial.data <- data
  final_res$nimput <- nmiss_initial
    
  if (model %in% c("all","aic","bic")) {
    final_res$tuning.parameters <- formulas
  }
  if (model %in% c("lasso", "ridge")) {
    final_res$tuning.parameters=list(lambda=lambdas)
  } else if (model == "elasticnet") {
    final_res$tuning.parameters=list(alpha=alphas, lambda=lambdas)
  }
  
  
  return(final_res)
}
