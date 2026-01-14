transport <- function(object, newdata, effect="ATE", boot.number=500, seed=NULL) {
  if (!inherits(object, c("gcbinary", "gctimes"))) {
    stop("object must be of class 'gcbinary' or 'gctimes'")
  }
  if (!is.null(object$newdata)) {stop("Cannot transport an already transported object")}
  
  fit <- object$calibration$fit 
  model <- object$model
  formula <- object$tuning.parameters
  all_terms <- attr(terms(object$formula), "term.labels")
  group <- object$group
  boot.type <- object$boot.type 
  
  
  newdata <- newdata[,which(colnames(newdata) %in% c(all_terms))]
  
  if (any(is.na(newdata))){
    nmiss <- nrow(newdata)
    newdata <- na.omit(newdata)
    nmiss <- nmiss - nrow(newdata)
    warning("Rows containing NA values have been removed from the dataset!")
  } else {
    nmiss <- 0
  }
  
  if(is.null(seed)) seed <- sample(1:10000,1)
  set.seed(seed)
  
  if (inherits(object, "gcbinary")) {
    outcome <- as.character(object$formula[[2]])
    if(model %in% c("lasso","ridge","elasticnet")) {
      
      stop("Not available for penalized methods.")
      
    } else { #glm
      beta.hat <- coef(fit)
      V.beta <- vcov(fit)
    }
    
    
    sim_betas <- MASS::mvrnorm(n = boot.number, mu = beta.hat, Sigma = V.beta)
    
    
    p0 <- c()
    p1 <- c()
    OR <- c()
    delta <- c()
    ratio <- c()
    
    for(b in 1:boot.number) {
      coef.mc <- sim_betas[b,]
      data.valid <- newdata
      
      if(effect == "ATE") {
        data.valid0 <- data.valid1 <- data.valid
      } else if(effect == "ATT") {
        data.valid0 <- data.valid1 <- data.valid[data.valid[,group] == 1,]
      } else { # ATU
        data.valid0 <- data.valid1 <- data.valid[data.valid[,group] == 0,]
      }
      
      data.valid0[,group] <- 0
      data.valid1[,group] <- 1
      
      if(model %in% c("lasso","ridge","elasticnet")) {
        .x0 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid0)[,-1]
        .x1 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid1)[,-1]
        .lp0 <- .x0 %*% coef.mc
        .lp1 <- .x1 %*% coef.mc
        .p0 <- mean(plogis(.lp0))
        .p1 <- mean(plogis(.lp1))
      } else {
        .X0 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid0)
        .X1 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid1)
        .lp0 <- .X0 %*% coef.mc
        .lp1 <- .X1 %*% coef.mc
        .p0 <- mean(plogis(.lp0))
        .p1 <- mean(plogis(.lp1))
      }
      
      .OR <- (.p1*(1-.p0))/(.p0*(1-.p1))
      .delta <- .p1 - .p0
      .ratio <- .p1 / .p0
      
      p0 <- c(p0, .p0)
      p1 <- c(p1, .p1)
      OR <- c(OR, .OR)
      delta <- c(delta, .delta)
      ratio <- c(ratio, .ratio)  
      }
    
    res <- list(calibration=object$calibration, 
                tuning.parameters=object$tuning.parameters, 
                data=object$data, 
                newdata=newdata,
                formula=formula, 
                model=model,
                cv=object$cv, 
                missing=nmiss,
                boot.number = boot.number,
                group=group,
                n = nrow(newdata),
                nevent = NA,
                adjusted.results = data.frame(p1 = p1, p0 = p0, delta = delta, ratio = ratio, OR = OR),
                effect=effect,
                call=match.call())
    class(res) <- "gcbinary"
    return(res)
  }
  
  
  if (inherits(object, "gctimes")) {
    
    
    
    
    if(model %in% c("lasso","ridge","elasticnet")) {
      stop("Not available for penalized methods.")
    }
    
    fit <- object$calibration$fit
    formula <- object$tuning.parameters
    group <- object$group
    boot.type <- object$boot.type
    
    beta.hat <- coef(fit)
    V.beta <- vcov(fit)
    
    sim_betas <- MASS::mvrnorm(n = boot.number, mu = beta.hat, Sigma = V.beta)
    
    AHR <- c()
    RMST0 <- c()
    RMST1 <- c()
    deltaRMST <- c()
    surv0 <- c()
    surv1 <- c()
    deltasurv <- c()
    
    for(b in 1:boot.number) {
      coef.mc <- sim_betas[b,]
      
      data.valid <- newdata
      if(effect == "ATE") {
        data.valid0 <- data.valid1 <- data.valid
      } else if(effect == "ATT") {
        data.valid0 <- data.valid1 <- data.valid[data.valid[,group] == 1,]
      } else { # ATU
        data.valid0 <- data.valid1 <- data.valid[data.valid[,group] == 0,]
      }
      
      data.valid0[,group] <- 0
      data.valid1[,group] <- 1
      
      .X0 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid0)[,-1]
      .X1 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid1)[,-1]
      
      .lp0 <- .X0 %*% coef.mc
      .lp1 <- .X1 %*% coef.mc
      
      H0.multi <- object$calibration$H0.multi
      T.multi <- object$calibration$time
      pro.time <- object$pro.time
      
      h0 <- diff(H0.multi)
      
      Si.0 <- exp(-exp(.lp0) %*% t(H0.multi))
      hi.0 <- exp(.lp0) %*% t(c(0,h0))
      h.mean.0 <- colSums(Si.0 * hi.0) / colSums(Si.0)
      H.mean.0 <- cumsum(h.mean.0)
      S.mean.0 <- exp(-H.mean.0)
      
      Si.1 <- exp(-exp(.lp1) %*% t(H0.multi))
      hi.1 <- exp(.lp1) %*% t(c(0,h0))
      h.mean.1 <- colSums(Si.1 * hi.1) / colSums(Si.1)
      H.mean.1 <- cumsum(h.mean.1)
      S.mean.1 <- exp(-H.mean.1)
      
      AHR <- c(AHR, sum(h.mean.1) / sum(h.mean.0))
      
      t.order <- order(T.multi)
      .S0_ord <- S.mean.0[t.order]
      .S1_ord <- S.mean.1[t.order]
      .T_ord <- T.multi[t.order]
      
      .t <- c(.T_ord[.T_ord <= pro.time], min(pro.time, max(.T_ord)))
      .s0 <- c(.S0_ord[.T_ord <= pro.time], .S0_ord[length(.S0_ord)])
      .s1 <- c(.S1_ord[.T_ord <= pro.time], .S1_ord[length(.S1_ord)])
      
      .RMST0 <- sum(diff(.t) * .s0[-length(.s0)])
      .RMST1 <- sum(diff(.t) * .s1[-length(.s1)])
      
      .surv0 <- .S0_ord[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
      .surv1 <- .S1_ord[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
      
      RMST0 <- c(RMST0, .RMST0)
      RMST1 <- c(RMST1, .RMST1)
      deltaRMST <- c(deltaRMST, .RMST1 - .RMST0)
      surv0 <- c(surv0, .surv0)
      surv1 <- c(surv1, .surv1)
      deltasurv <- c(deltasurv, .surv1 - .surv0)
    }
    
    res <- list(calibration = object$calibration,
                tuning.parameters = object$tuning.parameters,
                data = object$data,
                newdata = newdata,
                formula = formula,
                model = model,
                cv = object$cv,
                missing = nmiss,
                pro.time = pro.time,
                boot.number = boot.number,
                group = group,
                n = nrow(newdata),
                nevent = NA,
                adjusted.results = data.frame(AHR = AHR, RMST0 = RMST0, RMST1 = RMST1, deltaRMST = deltaRMST, s0 = surv0, s1 = surv1, delta = deltasurv),
                call = match.call())
    class(res) <- "gctimes"
    return(res)
  }
}