transport <- function(object, newdata, boot.number=500, seed=NULL) {
  if (!inherits(object, c("gcbinary", "gctimes", "gccount", "gccontinuous" ))) {
    stop("object must be of class 'gcbinary', 'gctimes', 'gccontinuous' or 'gccount'")
  }
  if (!is.null(object$newdata)) {stop("Cannot transport an already transported object")}
  
  fit <- object$calibration$fit 
  model <- object$model
  
  if(model %in% c("lasso","ridge","elasticnet")) {
    formula <- object$formula
  } else {
    formula <- object$tuning.parameters
  }
  
  all_terms <- attr(terms(formula), "term.labels")
  group <- object$group
  boot.type <- object$boot.type 
  
  
  
  all_vars <- all.vars(formula)
  if (inherits(object, "gctimes")) {
    all_vars <- all_vars[-1]
    all_vars <- all_vars[-1]
  } else {
    all_vars <- all_vars[-1]
  }
    
  all_vars <- setdiff(all_vars, object$group)
  
  missing_vars <- all_vars[!(all_vars %in% colnames(newdata))]
  
  if(length(missing_vars) > 0){
    stop(paste0("Some variables from the original model are not in newdata: ", paste(missing_vars, collapse=", ")))
  }
  
  
  if (any(is.na(newdata))){
    nmiss <- nrow(newdata)
    newdata <- na.omit(newdata)
    nmiss <- nmiss - nrow(newdata)
    warning("Rows containing NA values have been removed from the dataset!")
  } else {
    nmiss <- 0
  }
  
  if(is.null(seed)) {seed <- sample(1:10000,1)}
  set.seed(seed)
  
  if (inherits(object, "gcbinary")) {
    if(model %in% c("lasso","ridge","elasticnet")) {
      
      
      data.valid0 <- data.valid1 <- data.valid <- newdata
      
      data.valid0[,group] <- 0
      data.valid1[,group] <- 1
      
      .x0 <- model.matrix(update(formula, NULL ~ .), data.valid0)[,-1]
      .x1 <- model.matrix(update(formula, NULL ~ .), data.valid1)[,-1]
      
      .p0 <- mean(predict(fit, newx = .x0, type="response"))
      .p1 <- mean(predict(fit, newx = .x1, type="response"))
      
      .OR <- (.p1 * (1 - .p0)) / (.p0 * (1 - .p1))
      .delta <- .p1 - .p0
      .ratio <- .p1 / .p0
      
      res <- list(
        calibration = object$calibration,
        tuning.parameters = object$tuning.parameters,
        data = object$data,
        newdata = newdata,
        formula = formula,
        model = model,
        cv = object$cv,
        missing = nmiss,
        boot.number = 1,
        group = group,
        n = nrow(newdata) - nmiss,
        nevent = NA,
        adjusted.results = data.frame(p1 = .p1, p0 = .p0, delta = .delta, ratio = .ratio, OR = .OR),
        effect="ATE",
        call = match.call()
      )
      class(res) <- "gcbinary"
      return(res)
      
      
    } else { #glm
      beta.hat <- coef(fit)
      V.beta <- vcov(fit)

    sim_betas <- MASS::mvrnorm(n = boot.number, mu = beta.hat, Sigma = V.beta)
    
    
    p0 <- c()
    p1 <- c()
    OR <- c()
    delta <- c()
    ratio <- c()
    
    for(b in 1:boot.number) {
      coef.mc <- sim_betas[b,]
      
      data.valid0 <- data.valid1 <- data.valid <- newdata
      
      data.valid0[,group] <- 0
      data.valid1[,group] <- 1
      
        .X0 <- model.matrix(update(formula, NULL ~ .), data.valid0)
        .X1 <- model.matrix(update(formula, NULL ~ .), data.valid1)
        .lp0 <- .X0 %*% coef.mc
        .lp1 <- .X1 %*% coef.mc
        .p0 <- mean(plogis(.lp0))
        .p1 <- mean(plogis(.lp1))
      
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
                n = nrow(newdata) - nmiss,
                nevent = NA,
                adjusted.results = data.frame(p1 = p1, p0 = p0, delta = delta, ratio = ratio, OR = OR),
                effect="ATE",
                call=match.call())
    class(res) <- "gcbinary"
    return(res)
    }
  }
  
  
  if (inherits(object, "gctimes")) {
    H0.multi <- object$calibration$H0.multi
    T.multi <- object$calibration$time
    pro.time <- object$pro.time
    
    data.valid0 <- data.valid1 <- newdata
    data.valid0[, group] <- 0
    data.valid1[, group] <- 1
    
    if (model %in% c("all", "aic", "bic")) {
      .lp.0 <- predict(fit, newdata = data.valid0, type = "lp")
      .lp.1 <- predict(fit, newdata = data.valid1, type = "lp")
    } else {
      .x.valid0 <- model.matrix(update(formula, NULL ~ .), data.valid0)[,-1]
      .x.valid1 <- model.matrix(update(formula, NULL ~ .), data.valid1)[,-1]
      .lp.0 <- predict(fit, newx = .x.valid0)
      .lp.1 <- predict(fit, newx = .x.valid1)
    }
    
    lp.0 <- as.vector(.lp.0)
    lp.1 <- as.vector(.lp.1)
    
    h0 <- (H0.multi[2:length(T.multi)] - H0.multi[1:(length(T.multi)-1)])
    
    hi.0 <- exp(lp.0) * matrix(rep(h0, length(lp.0)), nrow = length(lp.0), byrow = TRUE)
    Si.0 <- exp(-exp(lp.0) * matrix(rep(H0.multi, length(lp.0)), nrow = length(lp.0), byrow = TRUE))
    hi.0 <- cbind(rep(0, length(lp.0)), hi.0)
    h.mean.0 <- apply(Si.0 * hi.0, FUN = "sum", MARGIN = 2) / apply(Si.0, FUN = "sum", MARGIN = 2)
    S.mean.0 <- exp(-cumsum(h.mean.0))
    
    hi.1 <- exp(lp.1) * matrix(rep(h0, length(lp.1)), nrow = length(lp.1), byrow = TRUE)
    Si.1 <- exp(-exp(lp.1) * matrix(rep(H0.multi, length(lp.1)), nrow = length(lp.1), byrow = TRUE))
    hi.1 <- cbind(rep(0, length(lp.1)), hi.1)
    h.mean.1 <- apply(Si.1 * hi.1, FUN = "sum", MARGIN = 2) / apply(Si.1, FUN = "sum", MARGIN = 2)
    S.mean.1 <- exp(-cumsum(h.mean.1))
    
    .AHR <- sum(h.mean.1) / sum(h.mean.0)
    
    .S.mean.0 <- S.mean.0[order(T.multi)]
    .S.mean.1 <- S.mean.1[order(T.multi)]
    .T.multi <- T.multi[order(T.multi)]
    
    .t <- c(.T.multi[.T.multi <= pro.time], min(pro.time, max(.T.multi)))
    .s0 <- c(.S.mean.0[.T.multi <= pro.time], .S.mean.0[length(.S.mean.0)])
    .s1 <- c(.S.mean.1[.T.multi <= pro.time], .S.mean.1[length(.S.mean.1)])
    
    .RMST0 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s0[1:(length(.s0) - 1)])
    .RMST1 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s1[1:(length(.s1) - 1)])
    
    .surv0 <- .S.mean.0[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
    .surv1 <- .S.mean.1[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
    
    res <- list(
      calibration = object$calibration,
      tuning.parameters = object$tuning.parameters,
      data = object$data,
      newdata = newdata,
      formula = formula,
      model = model,
      cv = object$cv,
      missing = nmiss,
      pro.time = pro.time,
      boot.number = 1,
      group = group,
      n = nrow(newdata) - nmiss,
      nevent = NA,
      adjusted.results = data.frame(AHR = .AHR, RMST0 = .RMST0, RMST1 = .RMST1, deltaRMST = .RMST1 - .RMST0, 
        s0 = .surv0, s1 = .surv1, delta = .surv1 - .surv0),
      effect="ATE",
      call = match.call()
    )
    class(res) <- "gctimes"
    return(res)
  }
  
  if (inherits(object, "gccount")) {
    if(model %in% c("lasso","ridge","elasticnet")) {
      data.valid0 <- data.valid1 <- data.valid <- newdata
      data.valid0[,group] <- 0
      data.valid1[,group] <- 1
      
      .x0 <- model.matrix(update(formula, NULL ~ .), data.valid0)[,-1]
      .x1 <- model.matrix(update(formula, NULL ~ .), data.valid1)[,-1]
      
      .c0 <- mean(predict(fit, newx = .x0, type = "response"))
      .c1 <- mean(predict(fit, newx = .x1, type = "response"))
      .delta <- .c1 - .c0
      .ratio <- .c1 / .c0
      
      res <- list(
        calibration = object$calibration,
        tuning.parameters = object$tuning.parameters,
        data = object$data,
        newdata = newdata,
        formula = formula,
        model = model,
        cv = object$cv,
        missing = nmiss,
        boot.number = 1, 
        group = group,
        n = nrow(newdata) - nmiss,
        adjusted.results = data.frame(c1 = .c1, c0 = .c0, delta = .delta, ratio = .ratio),
        effect="ATE",
        call = match.call()
      )
      class(res) <- "gccount" 
      return(res)
      
    } else { 
      beta.hat <- coef(fit)
      V.beta <- vcov(fit)
      sim_betas <- MASS::mvrnorm(n = boot.number, mu = beta.hat, Sigma = V.beta)
      
      c0 <- c(); c1 <- c(); delta <- c(); ratio <- c()
      
      for(b in 1:boot.number) {
        coef.mc <- sim_betas[b,]
        data.valid0 <- data.valid1 <- data.valid <- newdata
        data.valid0[,group] <- 0
        data.valid1[,group] <- 1
        
        .X0 <- model.matrix(update(formula, NULL ~ .), data.valid0)
        .X1 <- model.matrix(update(formula, NULL ~ .), data.valid1)
        
        .lp0 <- .X0 %*% coef.mc
        .lp1 <- .X1 %*% coef.mc
        
        .c0 <- mean(exp(.lp0)) 
        .c1 <- mean(exp(.lp1))
        
        .delta <- .c1 - .c0
        .ratio <- .c1 / .c0
        
        c0 <- c(c0, .c0); c1 <- c(c1, .c1)
        delta <- c(delta, .delta); ratio <- c(ratio, .ratio)
      }
      
      res <- list(
        calibration = object$calibration,
        tuning.parameters = object$tuning.parameters,
        data = object$data,
        newdata = newdata,
        formula = formula,
        model = model,
        cv = object$cv,
        missing = nmiss,
        boot.number = boot.number,
        group = group,
        n = nrow(newdata) - nmiss,
        adjusted.results = data.frame(c1 = c1, c0 = c0, delta = delta, ratio = ratio),
        effect = "ATE",
        call = match.call()
      )
      class(res) <- "gccount" 
      return(res)
    }
  }
  
  if (inherits(object, "gccontinuous")) {
    if (model %in% c("lasso", "ridge", "elasticnet")) {
      data.valid0 <- data.valid1 <- data.valid <- newdata 
      data.valid0[, group] <- 0 
      data.valid1[, group] <- 1 
      
      .x0 <- model.matrix(update(formula, NULL ~ .), data.valid0)[, -1] 
      .x1 <- model.matrix(update(formula, NULL ~ .), data.valid1)[, -1]
      
      .m0 <- mean(predict(fit, newx = .x0, type = "response"))
      .m1 <- mean(predict(fit, newx = .x1, type = "response")) 
      .delta <- .m1 - .m0
      .ratio <- .m1 / .m0 
      
      res <- list(
        calibration = object$calibration,
        tuning.parameters = object$tuning.parameters,
        data = object$data,
        newdata = newdata,
        formula = formula,
        model = model,
        cv = object$cv,
        missing = nmiss,
        boot.number = 1,
        group = group,
        n = nrow(newdata) - nmiss,
        adjusted.results = data.frame(m1 = .m1, m0 = .m0, delta = .delta, ratio = .ratio),
        effect="ATE",
        call = match.call()
      )
      class(res) <- "gccontinuous" 
      return(res) 
      
    } else { 
      beta.hat <- coef(fit) 
      V.beta <- vcov(fit) 
      sim_betas <- MASS::mvrnorm(n = boot.number, mu = beta.hat, Sigma = V.beta) 
      
      m0 <- c() 
      m1 <- c()
      delta <- c()
      ratio <- c() 
      
      for (b in 1:boot.number) {
        coef.mc <- sim_betas[b, ] 
        data.valid0 <- data.valid1 <- data.valid <- newdata 
        data.valid0[, group] <- 0 
        data.valid1[, group] <- 1 
        
        .X0 <- model.matrix(update(formula, NULL ~ .), data.valid0) 
        .X1 <- model.matrix(update(formula, NULL ~ .), data.valid1)
        
        .lp0 <- .X0 %*% coef.mc 
        .lp1 <- .X1 %*% coef.mc 
        
        .m0 <- mean(.lp0) 
        .m1 <- mean(.lp1) 
        
        .delta <- .m1 - .m0 
        .ratio <- .m1 / .m0 
        
        m0 <- c(m0, .m0)
        m1 <- c(m1, .m1) 
        delta <- c(delta, .delta)
        ratio <- c(ratio, .ratio) 
      }
      
      res <- list(
        calibration = object$calibration,
        tuning.parameters = object$tuning.parameters,
        data = object$data,
        newdata = newdata,
        formula = formula,
        model = model,
        cv = object$cv,
        missing = nmiss,
        boot.number = boot.number,
        group = group,
        n = nrow(newdata) - nmiss,
        adjusted.results = data.frame(m1 = m1, m0 = m0, delta = delta, ratio = ratio),
        effect = "ATE",
        call = match.call()
      )
      class(res) <- "gccontinuous"
      return(res)
    }
  }
}