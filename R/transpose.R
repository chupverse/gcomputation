transpose <- function(object, newdata, effect="ATE", boot.number=500, seed=NULL) {
  if (!inherits(object, c("gclogi", "gcsurv"))) {
    stop("object must be of class 'gclogi' or 'gcsurv'")
  }
  if (!is.null(object$newdata)) {stop("Cannot transpose an already transposed object")}
  
  fit <- object$calibration$fit 
  method <- object$method
  formula <- object$formula
  group <- object$group
  
  if (any(is.na(newdata))){
    nmiss <- nrow(newdata)
    newdata <- na.omit(newdata)
    nmiss <- nmiss - nrow(newdata)
    warning("Rows containing NA values have been removed from the dataset!")
  } else {
    nmiss <- 0
  }
  
  if(effect=="ATE"){ ttt <- which(newdata[,group] %in% c(0,1))
  } else if(effect=="ATT"){ ttt <- which(newdata[,group]==1)
  } else { ttt <- which(newdata[,group]==0) }
  newdata <- newdata[ttt,]
  
  if(is.null(seed)) seed <- sample(1:10000,1)
  set.seed(seed)
  
  if (inherits(object, "gclogi")) {
    outcome <- object$outcome
    boot.type <- object$boot.type 
    
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
    
    for(b in 1:boot.number) {
      id <- sample(1:nrow(newdata), size=nrow(newdata), replace=TRUE)
      data.learn = newdata[id,]
      
      if (boot.type == "bcv") {
        data.valid = newdata[-sort(unique(id)),]
        if (nrow(data.valid) == 0) {data.valid = newdata[id,]}
      } else {
        data.valid = newdata[id,]
      }
      
      if (effect == "ATE") {
        data.valid0 = data.valid1 = data.valid
      } else {
        if (effect == "ATT") {
          data.valid0 = data.valid1 = data.valid[data.valid[,group] == 1,]
        } else { #ATU
          data.valid0 = data.valid1 = data.valid[data.valid[,group] == 0,]
        }
      }
      data.valid0[,group] <- 0
      data.valid1[,group] <- 1
      

      fit.unadj <- glm(formula = as.formula(paste(outcome,"~",group)), data=data.learn, family="binomial")
      .p0.unadj = mean(predict(fit.unadj, newdata = data.valid0, type = "response"))
      .p1.unadj = mean(predict(fit.unadj, newdata = data.valid1, type = "response"))
      .OR.unadj = (.p1.unadj*(1-.p0.unadj))/(.p0.unadj*(1-.p1.unadj))
      .delta.unadj = .p1.unadj - .p0.unadj
      .ratio.unadj = .p1.unadj / .p0.unadj
      
      if(method %in% c("lasso","ridge","elasticnet")) {
        .x0 <- model.matrix(formula,data.valid0)[,-1]
        .x1 <- model.matrix(formula,data.valid1)[,-1]
        .p0 <- mean(predict(fit,newx=.x0,type="response"))
        .p1 <- mean(predict(fit,newx=.x1,type="response"))
      } else { 
        .p0 <- mean(predict(fit,newdata=data.valid0,type="response"))
        .p1 <- mean(predict(fit,newdata=data.valid1,type="response"))
      }
      .OR = (.p1*(1-.p0))/(.p0*(1-.p1))
      .delta = .p1 - .p0
      .ratio = .p1 / .p0
      
      p0 <- c(p0, .p0)
      p1 <- c(p1, .p1)
      OR <- c(OR, .OR)
      delta <- c(delta, .delta)
      ratio <- c(ratio, .ratio)
      p0.unadj <- c(p0.unadj, .p0.unadj)
      p1.unadj <- c(p1.unadj, .p1.unadj)
      OR.unadj <- c(OR.unadj, .OR.unadj)
      delta.unadj <- c(delta.unadj, .delta.unadj)
      ratio.unadj <- c(ratio.unadj, .ratio.unadj)
    }
    
    res <- list(calibration=object$calibration, 
                tuning.parameters=object$tuning.parameters, 
                data=object$data, 
                newdata=newdata,
                formula=formula, 
                method=method,
                cv=object$cv, 
                missing=nmiss,
                boot.number = boot.number,
                outcome=outcome,
                group=group,
                n = nrow(newdata),
                nevent = sum(newdata[,outcome]),
                p0=p0,
                p1=p1,
                delta=delta,
                ratio=ratio,
                OR=OR,
                p0.unadj=p0.unadj,
                p1.unadj=p1.unadj,
                delta.unadj=delta.unadj,
                ratio.unadj=ratio.unadj,
                OR.unadj=OR.unadj,
                effect=effect,
                call=match.call())
    class(res) <- "gclogi"
    return(res)
  }
  
  
  if (inherits(object, "gcsurv")) {
    fit <- object$calibration$model 
    H0.multi <- object$calibration$H0.multi 
    T.multi <- object$calibration$time
    times <- object$outcome$times
    failures <- object$outcome$failures
    pro.time <- object$pro.time
    boot.type <- object$boot.type 
    
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
    
    BCVerror <- 0 
    pro.time.extrapolate <- 0
    
    for (b in 1:boot.number) {
      id = sample(1:nrow(newdata), size = nrow(newdata), replace = TRUE)
      data.learn = newdata[id,]
      
      if (boot.type == "bcv") {data.valid = newdata[-sort(unique(id)),]
      if (nrow(data.valid) == 0) {data.valid = newdata[id,]}
      } else {
        data.valid = newdata[id,]
      }
      
      if (effect == "ATE") {
        data.valid0 = data.valid1 = data.valid
      } else {
        if (effect == "ATT") {
          data.valid0 = data.valid1 = data.valid[data.valid[,group] == 1,]
        } else { #ATU
          data.valid0 = data.valid1 = data.valid[data.valid[,group] == 0,]
        }
      }
      data.valid0[,group] <- 0
      data.valid1[,group] <- 1
      

      fit_unadj <- coxph(formula = as.formula(paste0("Surv(",times,",",failures,")~",group)), data=data.learn)
      .lp.0.unadj <- predict(fit_unadj, newdata = data.valid0, type="lp")
      .lp.1.unadj <- predict(fit_unadj, newdata = data.valid1, type="lp")
      .lp.unadj_learn <- predict(fit_unadj, newdata = data.learn, type="lp")
      .b.unadj <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.unadj_learn, centered = FALSE)
      hazard.unadj <- .b.unadj$cumulative_base_hazard
      fit_times.unadj <- .b.unadj$times
      H0.multi.unadj <- c(0, hazard.unadj[fit_times.unadj %in% sort(unique(data.learn[data.learn[,failures]==1,times]))] )
      T.multi.unadj <- c(0, fit_times.unadj[fit_times.unadj %in% sort(unique(data.learn[data.learn[,failures]==1,times]))] )
      lp.0.unadj <- as.vector(.lp.0.unadj)
      lp.1.unadj <- as.vector(.lp.1.unadj)
      h0.unadj <- (H0.multi.unadj[2:length(T.multi.unadj)] - H0.multi.unadj[1:(length(T.multi.unadj)-1)])
      
      Si.0.unadj <- exp(-exp(lp.0.unadj) * matrix(rep(H0.multi.unadj,length(lp.0.unadj)), nrow=length(lp.0.unadj), byrow=TRUE))
      hi.0.unadj <- exp(lp.0.unadj) * matrix(rep(h0.unadj,length(lp.0.unadj)), nrow=length(lp.0.unadj), byrow=TRUE)
      hi.0.unadj <- cbind(rep(0,length(lp.0.unadj)),hi.0.unadj)
      h.mean.0.unadj <- apply(Si.0.unadj * hi.0.unadj, FUN="sum", MARGIN=2) / apply(Si.0.unadj, FUN="sum", MARGIN=2)
      H.mean.0.unadj <- cumsum(h.mean.0.unadj)
      S.mean.0.unadj <- exp(-H.mean.0.unadj)
      
      Si.1.unadj <- exp(-exp(lp.1.unadj) * matrix(rep(H0.multi.unadj,length(lp.1.unadj)), nrow=length(lp.1.unadj), byrow=TRUE))
      hi.1.unadj <- exp(lp.1.unadj) * matrix(rep(h0.unadj,length(lp.1.unadj)), nrow=length(lp.1.unadj), byrow=TRUE)
      hi.1.unadj <- cbind(rep(0,length(lp.1.unadj)),hi.1.unadj)
      h.mean.1.unadj <- apply(Si.1.unadj * hi.1.unadj, FUN="sum", MARGIN=2) / apply(Si.1.unadj, FUN="sum", MARGIN=2)
      H.mean.1.unadj <- cumsum(h.mean.1.unadj)
      S.mean.1.unadj <- exp(-H.mean.1.unadj)
      
      if (max(T.multi.unadj) != max(fit_times.unadj)) {
        T.multi.unadj = c(T.multi.unadj, max(fit_times.unadj))
        H.mean.0.unadj = c(H.mean.0.unadj, max(H.mean.0.unadj))
        S.mean.0.unadj = c(S.mean.0.unadj, min(S.mean.0.unadj))
        H.mean.1.unadj = c(H.mean.1.unadj, max(H.mean.1.unadj))
        S.mean.1.unadj = c(S.mean.1.unadj, min(S.mean.1.unadj))
      }
      
      .AHR.unadj <- sum(h.mean.1.unadj) / sum(h.mean.0.unadj)
      .S.mean.0.unadj_ord <- S.mean.0.unadj[order(T.multi.unadj)]
      .S.mean.1.unadj_ord <- S.mean.1.unadj[order(T.multi.unadj)]
      .T.multi.unadj_ord <- T.multi.unadj[order(T.multi.unadj)]
      .t.unadj <- c(.T.multi.unadj_ord[.T.multi.unadj_ord <= pro.time], min(pro.time, max(.T.multi.unadj_ord)))
      .s0.unadj <- c(.S.mean.0.unadj_ord[.T.multi.unadj_ord <= pro.time], .S.mean.0.unadj_ord[length(.S.mean.0.unadj_ord)])
      .s1.unadj <- c(.S.mean.1.unadj_ord[.T.multi.unadj_ord <= pro.time], .S.mean.1.unadj_ord[length(.S.mean.1.unadj_ord)])
      .RMST0.unadj <- sum((.t.unadj[2:length(.t.unadj)] - .t.unadj[1:(length(.t.unadj) - 1)]) * .s0.unadj[1:(length(.s0.unadj) - 1)])
      .RMST1.unadj <- sum((.t.unadj[2:length(.t.unadj)] - .t.unadj[1:(length(.t.unadj) - 1)]) * .s1.unadj[1:(length(.s1.unadj) - 1)])
      .deltaRMST.unadj <- .RMST1.unadj - .RMST0.unadj
      .surv0.unadj <- .S.mean.0.unadj_ord[findInterval(pro.time, .T.multi.unadj_ord, rightmost.closed = TRUE)]
      .surv1.unadj <- .S.mean.1.unadj_ord[findInterval(pro.time, .T.multi.unadj_ord, rightmost.closed = TRUE)]
      .deltasurv.unadj <- .surv1.unadj - .surv0.unadj
      
      AHR.unadj <- c(AHR.unadj, .AHR.unadj)
      RMST0.unadj <- c(RMST0.unadj, .RMST0.unadj)
      RMST1.unadj <- c(RMST1.unadj, .RMST1.unadj)
      deltaRMST.unadj <- c(deltaRMST.unadj, .deltaRMST.unadj)
      surv0.unadj <- c(surv0.unadj, .surv0.unadj)
      surv1.unadj <- c(surv1.unadj, .surv1.unadj)
      deltasurv.unadj <- c(deltasurv.unadj, .deltasurv.unadj)
      
      if (max(data.learn[,times]) < pro.time) {
        pro.time.extrapolate = pro.time.extrapolate + 1
      }
      
      .lp.0 <- tryCatch({
        if (method %in% c("all", "aic", "bic")) { 
          predict(fit, newdata = data.valid0, type="lp")
        } else { 
          .x.valid0 <- model.matrix(formula, data.valid0)[,-1] 
          predict(fit, newx = .x.valid0)
        }
      }, error = function(e) {return(NULL)})
      if (is.null(.lp.0)) { BCVerror <- BCVerror + 1 ; next }
      
      .lp.1 <- tryCatch({
        if (method %in% c("all", "aic", "bic")) { 
          predict(fit, newdata = data.valid1, type="lp")
        } else { 
          .x.valid1 <- model.matrix(formula, data.valid1)[,-1] 
          predict(fit, newx = .x.valid1)
        }
      }, error = function(e) {return(NULL)})
      if (is.null(.lp.1)) { BCVerror <- BCVerror + 1 ; next }
      
      h0 <- (H0.multi[2:length(T.multi)] - H0.multi[1:(length(T.multi)-1)])
      
      lp.0 <- as.vector(.lp.0)
      lp.1 <- as.vector(.lp.1)
      
      Si.0 <- exp(-exp(lp.0) * matrix(rep(H0.multi,length(lp.0)), nrow=length(lp.0), byrow=TRUE))
      hi.0 <- exp(lp.0) * matrix(rep(h0,length(lp.0)), nrow=length(lp.0), byrow=TRUE)
      hi.0 <- cbind(rep(0,length(lp.0)),hi.0)
      h.mean.0 <- apply(Si.0 * hi.0, FUN="sum", MARGIN=2) / apply(Si.0, FUN="sum", MARGIN=2)
      H.mean.0 <- cumsum(h.mean.0)
      S.mean.0 <- exp(-H.mean.0)
      
      Si.1 <- exp(-exp(lp.1) * matrix(rep(H0.multi,length(lp.1)), nrow=length(lp.1), byrow=TRUE))
      hi.1 <- exp(lp.1) * matrix(rep(h0,length(lp.1)), nrow=length(lp.1), byrow=TRUE)
      hi.1 <- cbind(rep(0,length(lp.1)),hi.1)
      h.mean.1 <- apply(Si.1 * hi.1, FUN="sum", MARGIN=2) / apply(Si.1, FUN="sum", MARGIN=2)
      H.mean.1 <- cumsum(h.mean.1)
      S.mean.1 <- exp(-H.mean.1)
      
      if (max(T.multi) != max(object$calibration$time)) {
        T.multi = c(T.multi, max(object$calibration$time))
        H.mean.0 = c(H.mean.0, max(H.mean.0))
        S.mean.0 = c(S.mean.0, min(S.mean.0))
        H.mean.1 = c(H.mean.1, max(H.mean.1))
        S.mean.1 = c(S.mean.1, min(S.mean.1))
      }
      
      AHR <- c(AHR, sum(h.mean.1) / sum(h.mean.0) )
      
      .S.mean.0_ord <- S.mean.0[order(T.multi)]
      .S.mean.1_ord <- S.mean.1[order(T.multi)]
      .T.multi_ord <- T.multi[order(T.multi)]
      
      .t <- c(.T.multi_ord[.T.multi_ord <= pro.time], min(pro.time, max(.T.multi_ord)))
      .s0 <- c(.S.mean.0_ord[.T.multi_ord <= pro.time], .S.mean.0_ord[length(.S.mean.0_ord)])
      .s1 <- c(.S.mean.1_ord[.T.multi_ord <= pro.time], .S.mean.1_ord[length(.S.mean.1_ord)])
      
      .RMST0 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s0[1:(length(.s0) - 1)])
      .RMST1 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s1[1:(length(.s1) - 1)])
      .surv1 <- .S.mean.1_ord[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
      .surv0 <- .S.mean.0_ord[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
      
      surv1 <- c(surv1,.surv1)
      surv0 <- c(surv0,.surv0)
      deltasurv <- c(deltasurv, .surv1 - .surv0)
      RMST0 <- c(RMST0, .RMST0)
      RMST1 <- c(RMST1, .RMST1)
      deltaRMST <- c(deltaRMST, .RMST1 - .RMST0)
    }
    
    if (pro.time.extrapolate > 1) {warning(paste0("In at least one boostrap sample the \"pro.time\" was higher than the maximum follow-up time (survival was extrapolated in ",pro.time.extrapolate," bootstrap samples). It is advised to pick a lower value for \"pro.time\""))}
    if (BCVerror > 1) {warning(paste0("Skipped ",BCVerror," bootstrap iterations due to the validation dataset containing factors not in the train dataset. Either use type=\"boot\" instead of \"bcv\" or remove factors with rare modalities."))}
    
    res <- list(calibration=object$calibration,
                tuning.parameters=object$tuning.parameters, 
                data=object$data, 
                newdata=newdata,
                formula=formula, 
                method=method,
                cv=object$cv, 
                missing=nmiss,
                pro.time=pro.time,
                boot.number = boot.number,
                outcome=object$outcome,
                group=group,
                n = nrow(newdata),
                nevent = sum(newdata[,failures]),
                AHR = AHR,
                RMST0 = RMST0,
                RMST1 = RMST1,
                deltaRMST = deltaRMST,
                surv0 = surv0,
                surv1 = surv1,
                deltasurv = deltasurv,
                AHR.unadj = AHR.unadj,
                RMST0.unadj = RMST0.unadj,
                RMST1.unadj = RMST1.unadj,
                deltaRMST.unadj = deltaRMST.unadj,
                surv0.unadj = surv0.unadj,
                surv1.unadj = surv1.unadj,
                deltasurv.unadj = deltasurv.unadj,
                call = match.call())
    class(res) <- "gcsurv"
    return(res)
  }
}