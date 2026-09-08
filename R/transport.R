transport <- function(object, newdata, estim_var, nboot=1000, n.sim=500, seed=NULL) {
  if (!inherits(object, c("gcbinary", "gctimes", "gccount", "gccontinuous" ))) {
    stop("object must be of class 'gcbinary', 'gctimes', 'gccontinuous' or 'gccount'")
  }
  if (!is.null(object$newdata)) {stop("Cannot transport an already transported object")}
  
  if("initial.data" %in% attributes(object)$names){stop("Cannot transport when multiple imputations has been used, try relaunch gcomputation without it.")}
  
  if(!(estim_var %in% c("monte_carlo", "point_estimate", "m_estimation", "bootstrap"))){
    stop("estim_ var parameter needs to be one of: monte_carlo, point_estimate, m_estimation, bootstrap")
  }
  
  fit <- object$qmodel.fit
  model <- object$model
  
  if(estim_var == "m_estimation"){
    if(!(model %in% c("all", "aic", "bic"))){
      stop("M-estimation is only defined for parametric regression models: model = all, aic, bic")
    }
  }
  
  if(estim_var == "monte_carlo"){
    if(!(model %in% c("all", "aic", "bic"))){
      stop("Monte Carlo simulation is only defined for parametric regression models: model = all, aic, bic")
    }
  }
  
  if(estim_var == "point_estimate"){
    if(!((model %in% c("lasso", "ridge", "elasticnet", "")) || inherits(object, "gctimes"))){
      stop("Point estimate is only defined for penalized models (lasso, ridge, elasticnet) and survival models.")
    }
  }
  
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
    warning("Rows containing NA values have been removed from the targeted dataset!")
  } else {
    nmiss <- 0
  }
  
  if(!is.null(seed)) {set.seed(seed)}
  set.seed(seed)
  
  if(estim_var %in% c("monte_carlo", "point_estimate")){
    
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
          qmodel.fit = object$qmodel.fit,
          predictions = object$predictions,
          tuning.parameters = object$tuning.parameters,
          data = object$data,
          newdata = newdata,
          formula = formula,
          model = model,
          cv = object$cv,
          missing = nmiss,
          n.sim = 1,
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
        
        sim_betas <- MASS::mvrnorm(n = n.sim, mu = beta.hat, Sigma = V.beta)
        
        
        p0 <- c()
        p1 <- c()
        OR <- c()
        delta <- c()
        ratio <- c()
        
        for(b in 1:n.sim) {
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
        
        res <- list(qmodel.fit = object$qmodel.fit,
                    predictions = object$predictions,
                    tuning.parameters=object$tuning.parameters, 
                    data=object$data, 
                    newdata=newdata,
                    formula=formula, 
                    model=model,
                    cv=object$cv, 
                    missing=nmiss,
                    n.sim = n.sim,
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
        qmodel.fit = object$qmodel.fit,
        calibration = object$calibration,
        tuning.parameters = object$tuning.parameters,
        data = object$data,
        newdata = newdata,
        formula = formula,
        model = model,
        cv = object$cv,
        missing = nmiss,
        pro.time = pro.time,
        n.sim = 1,
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
          n.sim = 1, 
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
        sim_betas <- MASS::mvrnorm(n = n.sim, mu = beta.hat, Sigma = V.beta)
        
        c0 <- c(); c1 <- c(); delta <- c(); ratio <- c()
        
        for(b in 1:n.sim) {
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
          n.sim = n.sim,
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
          qmodel.fit = object$qmodel.fit,
          predictions = object$predictions,
          tuning.parameters = object$tuning.parameters,
          data = object$data,
          newdata = newdata,
          formula = formula,
          model = model,
          cv = object$cv,
          missing = nmiss,
          n.sim = 1,
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
        sim_betas <- MASS::mvrnorm(n = n.sim, mu = beta.hat, Sigma = V.beta) 
        
        m0 <- c() 
        m1 <- c()
        delta <- c()
        ratio <- c() 
        
        for (b in 1:n.sim) {
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
          qmodel.fit = object$qmodel.fit,
          predictions = object$predictions,
          tuning.parameters = object$tuning.parameters,
          data = object$data,
          newdata = newdata,
          formula = formula,
          model = model,
          cv = object$cv,
          missing = nmiss,
          n.sim = n.sim,
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
  
  if(estim_var == "m_estimation"){
    
    if(inherits(object, "gctimes")){
      stop("M-estimation is not available for gc_times objects for now.")
    }
    
    list_mest <- create_mestim_obj(gc = object, data_target = newdata)
    res_mest <- Mestimation_process(list_mest)
    
    
    
    if(inherits(object, "gcbinary")){
      
      adj_results <- data.frame(p1 = res_mest["tau1"], 
                                p0 = res_mest["tau0"], 
                                delta = res_mest["delta"], 
                                ratio = res_mest["ratio"],
                                OR = res_mest["OR"],
                                sd_p1 = res_mest["sd_tau1"],
                                sd_p0 = res_mest["sd_tau0"],
                                sd_delta = res_mest["sd_delta"],
                                sd_ratio = res_mest["sd_ratio"],
                                sd_OR = res_mest["sd_OR"])
      
    }
    
    if(inherits(object, "gccontinuous")){
      
      adj_results <- data.frame(m1 = res_mest["tau1"], 
                                m0 = res_mest["tau0"], 
                                delta = res_mest["delta"], 
                                ratio = res_mest["ratio"],
                                sd_m1 = res_mest["sd_tau1"],
                                sd_m0 = res_mest["sd_tau0"],
                                sd_delta = res_mest["sd_delta"],
                                sd_ratio = res_mest["sd_ratio"])
      
    }
    
    if(inherits(object, "gccount")){
      
      adj_results <- data.frame(c1 = res_mest["tau1"], 
                                c0 = res_mest["tau0"], 
                                delta = res_mest["delta"], 
                                ratio = res_mest["ratio"],
                                sd_c1 = res_mest["sd_tau1"],
                                sd_c0 = res_mest["sd_tau0"],
                                sd_delta = res_mest["sd_delta"],
                                sd_ratio = res_mest["sd_ratio"])
      
    }
    
    .mm <- attr(res_mest, "model.matrix")
    new_qmodel.fit <- list()
    new_qmodel.fit$coefficients <- setNames(res_mest[1:length(colnames(.mm))],
                                            nm = colnames(.mm))
    
    res <- list(
      qmodel.fit = object$qmodel.fit,
      new_qmodel.fit = new_qmodel.fit,
      predictions = NA,
      tuning.parameters = object$tuning.parameters,
      data = object$data,
      newdata = newdata,
      formula = formula,
      model = model,
      cv = NULL,
      missing = nmiss,
      missing_origin = list_mest$nmiss_origin,
      n.sim = NULL,
      group = group,
      n = nrow(newdata) - nmiss,
      nevent = NA,
      adjusted.results = adj_results,
      effect = "ATE",
      call = match.call()
    )
    
    class(res) <- class(object)
    attr(res, "estim_var") <- "m_estimation"
    
    return(res)
    
  }
  
  if(estim_var == "bootstrap"){
    
    if(inherits(object, "gctimes")){
      stop("bootstrap for object gctimes not implemented.")
    }
    
    if(model %in% c("aic", "bic")){
      stop("bootstrap is not yet implemented for models selected using AIC or BIC.")
    }
    
    form <- object$formula
    data_form <- object$data %>%
      dplyr::select(all.vars(form))  
    
    if (any(is.na(data_form))){
      
      initial_data_omit <- na.omit(data_form)
      nmiss_origin <- nrow(data_form) - nrow(initial_data_omit)
      
      data_origin <- initial_data_omit
      
      warning("Rows containing NA values in the original dataset have been removed!")
      
    } else {
      
      data_origin <- data_form
      nmiss_origin <- 0
      
    }
    
    if(model %in% c("lasso", "ridge", "elasticnet")){
      
      list_res <- replicate(nboot,
                            penalized_transport_forboot(gc = object, 
                                                        data_ini = data_origin,
                                                        data_target = newdata,
                                                        boot = T),
                            simplify = F)
      
      res_boot <- dplyr::bind_rows(list_res)
      
      
    }
    
    if(model == "all"){
      
      list_res <- replicate(nboot,
                            parametric_transport_forboot(gc = object, 
                                                         data_ini = data_origin,
                                                         data_target = newdata,
                                                         boot = T),
                            simplify = F)
      
      res_boot <- dplyr::bind_rows(list_res)
      
      
    }
    
    if(inherits(object, "gcbinary")){
      
      ratio <- res_boot[["tau1"]]/res_boot[["tau0"]]
      OR <- (res_boot[["tau1"]]*(1 - res_boot[["tau0"]]))/(res_boot[["tau0"]]*(1 - res_boot[["tau1"]]))
      
      adj_results <- data.frame(p1 = res_boot[["tau1"]], 
                                p0 = res_boot[["tau0"]], 
                                delta = res_boot[["delta"]], 
                                ratio = ratio,
                                OR = OR)
      
    }
    
    if(inherits(object, "gccontinuous")){
      
      ratio <- res_boot[["tau1"]]/res_boot[["tau0"]]
      
      adj_results <- data.frame(m1 = res_boot[["tau1"]], 
                                m0 = res_boot[["tau0"]], 
                                delta = res_boot[["delta"]], 
                                ratio = ratio)
      
      
    }
    
    if(inherits(object, "gccount")){
      
      ratio <- res_boot[["tau1"]]/res_boot[["tau0"]]
      
      adj_results <- data.frame(c1 = res_boot[["tau1"]], 
                                c0 = res_boot[["tau0"]], 
                                delta = res_boot[["delta"]], 
                                ratio = ratio)
      
    }
    
    
    
    res <- list(
      qmodel.fit = object$qmodel.fit,
      predictions = NA,
      tuning.parameters = object$tuning.parameters,
      data = object$data,
      newdata = newdata,
      formula = formula,
      model = model,
      cv = object$cv,
      missing = nmiss,
      missing_origin = nmiss_origin,
      n.sim = NULL,
      nboot = nboot,
      group = group,
      n = nrow(newdata) - nmiss,
      nevent = NA,
      adjusted.results = adj_results,
      effect = "ATE",
      call = match.call()
    )
    
    class(res) <- class(object)
    
    return(res)
    
    
  }
  
  
}

create_mestim_obj <- function(gc, data_target){
  
  form <- gc$tuning.parameters
  
  data_form <- gc$data %>%
    dplyr::select(all.vars(form))  
  
  if (any(is.na(data_form))){
    
    initial_data_omit <- na.omit(data_form)
    nmiss_origin <- nrow(data_form) - nrow(initial_data_omit)
    
    data_origin <- initial_data_omit
    
    warning("Rows containing NA values in the original dataset have been removed!")
    
  } else {
    
    data_origin <- data_form
    nmiss_origin <- 0
    
  }
  
  grp_var <- gc$group
  family <- gc$qmodel.fit$family$family
  
  data_origin$S <- 1
  data_target$S <- 0
  
  outcome_var <- all.vars(form)[1]
  data_target[ , outcome_var] <- NA 
  
  data_all <- bind_rows(data_origin, data_target) %>%
    dplyr::select(intersect(names(data_origin), names(data_target)))
  
  n <- sum(data_all$S) #number of individuals in the original data set
  m <- nrow(data_all) - n #number of individuals in the targeted data set
  
  root_start <- c(gc$qmodel.fit$coefficients,
                  mean(gc$adjusted.results[,2]),
                  mean(gc$adjusted.results[,1]),
                  mean(gc$adjusted.results$delta),
                  mean(gc$adjusted.results$ratio))
  
  
  if(family == "binomial"){
    root_start <- c(root_start, mean(gc$adjusted.results$OR))
  }
  
  out <- list(data = data_all,
              n = n,
              m = m,
              nmiss_origin = nmiss_origin,
              grp_var = grp_var,
              formula = form,
              family = family,
              root_start = root_start,
              qmodel_fit = gc$qmodel.fit)
  
  class(out) <- c(paste0("mestim_", family), "mestim_list")
  
  return(out)
}

Mestimation_process.default <- function(mestim_list, ...){
  stop("Family not supported for M-estimation: ", mestim_list$family)
}

Mestimation_process.mestim_binomial <- function(mestim_list, ...){
  link_fun <- function(eta) plogis(eta)
  .mestim_process_core(mestim_list, 
                       link_fun = link_fun,
                       estimands = list(
                         delta = function(p0,p1){p1 - p0},
                         ratio = function(p0, p1){p1/p0},
                         OR = function(p0, p1){(p1 * (1 - p0)) / (p0 * (1 - p1))}
                       ))
}

Mestimation_process.mestim_poisson <- function(mestim_list, ...){
  link_fun <- function(eta) exp(eta)
  .mestim_process_core(mestim_list, 
                       link_fun = link_fun,
                       estimands = list(
                         delta = function(c0,c1){c1 - c0},
                         ratio = function(c0, c1){c1/c0}
                       ))
}

Mestimation_process.mestim_gaussian <- function(mestim_list, ...){
  link_fun <- function(eta) eta
  .mestim_process_core(mestim_list, 
                       link_fun = link_fun,
                       estimands = list(
                         delta = function(m0,m1){m1 - m0},
                         ratio = function(m0, m1){m1/m0}
                       ))
}

.mestim_process_core <- function(mestim_list, link_fun, estimands){
  
  n_estimands    <- length(estimands)
  estimands_names <- names(estimands)
  
  fun_mestimate <- function(data, m, n, formula, model_terms, group_var, family){
    
    vars_f <- all.vars(formula)
    
    if(is.factor(data[, vars_f[1]])){
      
      Y <- as.numeric(as.character(data[, vars_f[1]]))
      
    }else{
      
      Y <- as.numeric(data[, vars_f[1]])
      
    }
    
    Y[is.na(Y)] <- 0 #to avoid computational issues due to NA in m_estimate function (see (*) below)
    
    A  <- as.numeric(as.character(data[ ,group_var])) 
    S  <- data$S
    
    data0 <- data1 <- data
    data0[ ,group_var] <- 0
    data1[ ,group_var] <- 1
    
    mm  <- model.matrix(model_terms, data = data)
    mm0 <- model.matrix(model_terms, data = data0)
    mm1 <- model.matrix(model_terms, data = data1)
    
    
    
    function(theta){
      
      p <- length(theta) - (2 + n_estimands)  # +2 due to tau0 and tau1 which appear in all cases
      beta <- theta[1:p]
      
      lp <- beta %*% t(mm)
      lp0 <- beta %*% t(mm0)
      lp1 <- beta %*% t(mm1) 
      
      out <- c(
        (S*(Y - link_fun(lp))) %*% mm, # (*) NA <- 0 useful here, Y with NA are the ones with S = 0 and would not be consider here.
        ((m + n)/m)*(1-S)*link_fun(lp0) - theta[p+1],
        ((m + n)/m)*(1-S)*link_fun(lp1) - theta[p+2]
      )
      
      
      
      eq_estimands <- vapply(seq_len(n_estimands), function(i){
        estimands[[i]](theta[p+1], theta[p+2]) - theta[p+2+i]
      }, 1)
      
      out <- c(out, eq_estimands)
      
      return(out)
      
    }
    
  }
  
  results <- geex::m_estimate(
    estFUN = fun_mestimate, 
    data   = mestim_list$data,
    outer_args = list(n = mestim_list$n, 
                      m = mestim_list$m, 
                      formula = mestim_list$form,
                      model_terms = delete.response(terms(mestim_list$qmodel_fit)), 
                      group_var = mestim_list$grp_var,
                      family = mestim_list$family),
    root_control = geex::setup_root_control(start = mestim_list$root_start))
  
  estim <- geex::coef(results)
  sd_estim <- sqrt(diag(geex::vcov(results)))
  
  mm <- model.matrix(object = delete.response(terms(mestim_list$formula)), data = mestim_list$data) #only to get correct variables names
  
  base_names <- c(colnames(mm), "tau0", "tau1")
  all_names  <- c(base_names, estimands_names)
  
  out <- setNames(
    c(estim, sd_estim),
    c(all_names, paste0("sd_", all_names))
  )
  
  attr(out, "model.matrix") <- mm
  
  return(out)
  
}

Mestimation_process <- function(mestim_list, ...){
  UseMethod("Mestimation_process")
}

penalized_transport_forboot <- function(gc, data_ini, data_target, boot = T){
  
  grp <- gc$group
  form <- gc$formula
  model <- gc$model
  
  if(isTRUE(boot)){
    
    data_train <- dplyr::slice_sample(.data = data_ini, 
                                      n = nrow(data_ini), 
                                      replace = T)
    
    data_test <- dplyr::slice_sample(.data = data_target, 
                                     n = nrow(data_target), 
                                     replace = T)
  }else{
    
    data_train <- data_ini
    data_test <- data_target
    
  }
  
  data_test0 <- data_test1 <- data_test
  data_test0[,grp] <- 0
  data_test1[,grp] <- 1
  
  mm <- model.matrix(object = form, data = data_train)
  mm0 <- model.matrix(object = form, data = data_test0)
  mm1 <- model.matrix(object = form, data = data_test1)
  
  if("(Intercept)" %in% colnames(mm)){
    
    mm <- mm[,-1]
    mm0 <- mm0[,-1]
    mm1 <- mm1[,-1]
    
    intercept_bool <- T
    
  }else{
    
    intercept_bool <- F
    
  }
  
  y <- data_train[,all.vars(form)[1]]
  
  if(model == "ridge"){
    alpha <- 0
  }else{
    if(model == "lasso"){
      alpha <- 1
    }else{
      alpha <- gc$tuning.parameters$alpha
    }
  }
  
  if(inherits(gc, "gcbinary")){
    fam <- "binomial"
  }else{
    if(inherits(gc, "gccontinuous")){
      fam <- "gaussian"
    }else{
      if(inherits(gc, "gccount")){
        fam <- "poisson"
      }
    }
  }
  
  penalty_fact <- rep(1, ncol(mm))
  penalty_fact[which(colnames(mm) == grp)] <- 0 #to set no penalty on group variable
  
  res_glmnet <- glmnet::glmnet(x = mm, y = y,family = fam, 
                               alpha = alpha, 
                               lambda = gc$tunning.parameters$lambda,
                               penalty.factor = penalty_fact, 
                               intercept = intercept_bool)
  
  
  res0 <- predict(res_glmnet, newx = mm0, type = "response")
  res1 <- predict(res_glmnet, newx = mm1, type = "response")
  
  tau0 <- mean(res0)
  tau1 <- mean(res1)
  ate <- tau1-tau0
  
  out <- c("tau0" = tau0,
           "tau1" = tau1,
           "delta" = ate)
  
  return(out)
}

parametric_transport_forboot <- function(gc, data_ini, data_target, boot = T){
  
  grp <- gc$group
  form <- gc$formula
  model <- gc$model
  
  if(isTRUE(boot)){
    
    data_train <- dplyr::slice_sample(.data = data_ini, 
                                      n = nrow(data_ini), 
                                      replace = T)
    
    data_test <- dplyr::slice_sample(.data = data_target, 
                                     n = nrow(data_target), 
                                     replace = T)
  }else{
    
    data_train <- data_ini
    data_test <- data_target
    
  }
  
  data_test0 <- data_test1 <- data_test
  data_test0[,grp] <- 0
  data_test1[,grp] <- 1
  
  
  if(inherits(gc, "gcbinary")){
    fam <- "binomial"
  }else{
    if(inherits(gc, "gccontinuous")){
      fam <- "gaussian"
    }else{
      if(inherits(gc, "gccount")){
        fam <- "poisson"
      }
    }
  }
  
  res_glm <- glm(formula = form,
                 family = fam, 
                 data = data_train)
  
  
  res0 <- predict(res_glm, newdata = data_test0, type = "response")
  res1 <- predict(res_glm, newdata = data_test1, type = "response")
  
  tau0 <- mean(res0)
  tau1 <- mean(res1)
  ate <- tau1-tau0
  
  out <- c("tau0" = tau0,
           "tau1" = tau1,
           "delta" = ate)
  
  return(out)
}

