miboot_gc_survival <- function(formula, data, group, pro.time, ..., m = 5,
                               effect = "ATE", gc.method, param.tune = NULL, cv = 10,
                               boot.type = "bcv", boot.number = 500,
                               boot.tune = FALSE, progress = TRUE, seed = NULL) {

  cl <- match.call()
  
  if (missing(pro.time)) {
    pro.time <- NULL
  }
  
  
  
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
  surv0 <- c()
  surv1 <- c()
  deltasurv <- c()
  AHR.unadj <- c()
  RMST0.unadj <- c()
  RMST1.unadj <- c()
  deltaRMST.unadj <- c()
  surv0.unadj <- c()
  surv1.unadj <- c()
  deltasurv.unadj <- c()
  
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
    
    gc_res <- gc_survival(formula = formula, data = current_imputed_data, group = group,
                          pro.time = pro.time, effect = effect, method = gc.method, param.tune = param.tune,
                          cv = cv, boot.type = boot.type, boot.number = boot.number,
                          boot.tune = boot.tune, progress = FALSE, seed = gc_seed)
    
    datas[[i]] <- current_imputed_data
    calibrations[[i]] <- gc_res$calibration
    
    AHR <- c(AHR, gc_res$AHR)
    RMST0 <- c(RMST0, gc_res$RMST0)
    RMST1 <- c(RMST1, gc_res$RMST1)
    deltaRMST <- c(deltaRMST, gc_res$deltaRMST)
    surv0 <- c(surv0, gc_res$surv0)
    surv1 <- c(surv1, gc_res$surv1)
    deltasurv <- c(deltasurv, gc_res$deltasurv)
    AHR.unadj <- c(AHR.unadj, gc_res$AHR.unadj)
    RMST0.unadj <- c(RMST0.unadj, gc_res$RMST0.unadj)
    RMST1.unadj <- c(RMST1.unadj, gc_res$RMST1.unadj)
    deltaRMST.unadj <- c(deltaRMST.unadj, gc_res$deltaRMST.unadj)
    surv0.unadj <- c(surv0.unadj, gc_res$surv0.unadj)
    surv1.unadj <- c(surv1.unadj, gc_res$surv1.unadj)
    deltasurv.unadj <- c(deltasurv.unadj, gc_res$deltasurv.unadj)
    
    if (gc.method %in% c("lasso", "ridge")) {
      if (!is.null(gc_res$tuning.parameters$lambda)) {
        lambdas <- c(lambdas, gc_res$tuning.parameters$lambda)
      }
    } else if (gc.method == "elasticnet") {
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
  
  final_res$AHR <- AHR
  final_res$RMST0 <- RMST0
  final_res$RMST1 <- RMST1
  final_res$deltaRMST <- deltaRMST
  final_res$surv0 <- surv0
  final_res$surv1 <- surv1
  final_res$deltasurv <- deltasurv
  final_res$AHR.unadj <- AHR.unadj
  final_res$RMST0.unadj <- RMST0.unadj
  final_res$RMST1.unadj <- RMST1.unadj
  final_res$deltaRMST.unadj <- deltaRMST.unadj
  final_res$surv0.unadj <- surv0.unadj
  final_res$surv1.unadj <- surv1.unadj
  final_res$deltasurv.unadj <- deltasurv.unadj
  
  final_res$boot.number <- m * boot.number
  final_res$seed <- seed
  final_res$call <- cl
  final_res$m <- m
  final_res$initial.data <- data
  
  if (gc.method %in% c("lasso", "ridge")) {
    final_res$tuning.parameters <- list(lambda = lambdas)
  } else if (gc.method == "elasticnet") {
    final_res$tuning.parameters <- list(alpha = alphas, lambda = lambdas)
  }
  
  return(final_res)
}
