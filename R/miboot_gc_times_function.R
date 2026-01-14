.miboot_gc_times <- function(formula, data, group, pro.time=NULL,
                               effect = "ATE", model, param.tune = NULL, cv = 10,
                               boot.type = "bcv", boot.number = 500,
                               boot.tune = FALSE, progress = TRUE, seed = NULL, m = 5, ...) {

  cl <- match.call()
  

  
  
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  mice_args <- list(...)
  mice_args$data <- data
  mice_args$m <- m
  mice_args$printFlag <- FALSE
  
  imp_object <- do.call(mice::mice, mice_args)
  
  AHR <- c()
  RMST0 <- c()
  RMST1 <- c()
  deltaRMST <- c()
  s0 <- c()
  s1 <- c()
  delta <- c()
  AHR.unadj <- c()
  RMST0.unadj <- c()
  RMST1.unadj <- c()
  deltaRMST.unadj <- c()
  s0.unadj <- c()
  s1.unadj <- c()
  delta.unadj <- c()
  
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
    
    gc_res <- .gc_times(formula = formula, data = current_imputed_data, group = group,
                          pro.time = pro.time, effect = effect, model = model, param.tune = param.tune,
                          cv = cv, boot.type = boot.type, boot.number = boot.number,
                          boot.tune = boot.tune, progress = FALSE, seed = gc_seed)
    
    datas[[i]] <- current_imputed_data
    calibrations[[i]] <- gc_res$calibration
    
    AHR <- c(AHR, gc_res$adjusted.results$AHR)
    RMST0 <- c(RMST0, gc_res$adjusted.results$RMST0)
    RMST1 <- c(RMST1, gc_res$adjusted.results$RMST1)
    deltaRMST <- c(deltaRMST, gc_res$adjusted.results$deltaRMST)
    s0 <- c(s0, gc_res$adjusted.results$s0)
    s1 <- c(s1, gc_res$adjusted.results$s1)
    delta <- c(delta, gc_res$adjusted.results$delta)
    AHR.unadj <- c(AHR.unadj, gc_res$unadjusted.results$AHR)
    RMST0.unadj <- c(RMST0.unadj, gc_res$unadjusted.results$RMST0)
    RMST1.unadj <- c(RMST1.unadj, gc_res$unadjusted.results$RMST1.unadj)
    deltaRMST.unadj <- c(deltaRMST.unadj, gc_res$unadjusted.results$deltaRMST)
    s0.unadj <- c(s0.unadj, gc_res$unadjusted.results$s0)
    s1.unadj <- c(s1.unadj, gc_res$unadjusted.results$s1)
    delta.unadj <- c(delta.unadj, gc_res$unadjusted.results$delta)
    
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
  
  final_res$adjusted.results <- data.frame(AHR = AHR, RMST0 = RMST0, RMST1 = RMST1, deltaRMST = deltaRMST, s0 = s0, s1 = s1, delta = delta)
  final_res$unadjusted.results <- data.frame(AHR = AHR.unadj, RMST0 = RMST0.unadj, RMST1 = RMST1.unadj, deltaRMST = deltaRMST.unadj, s0 = s0.unadj, s1 = s1.unadj, delta = delta.unadj)
  
  
  final_res$boot.number <- m * boot.number
  final_res$seed <- seed
  final_res$call <- cl
  final_res$m <- m
  final_res$initial.data <- data
  
  if (model %in% c("all","aic","bic")) {
    final_res$tuning.parameters <- formulas
  }
  if (model %in% c("lasso", "ridge")) {
    final_res$tuning.parameters <- list(lambda = lambdas)
  } else if (model == "elasticnet") {
    final_res$tuning.parameters <- list(alpha = alphas, lambda = lambdas)
  }
  
  return(final_res)
}
