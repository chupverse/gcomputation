transport <- function(object, newdata, effect="ATE", boot.number=500, seed=NULL) {
  if (!inherits(object, c("gcbinary", "gctimes"))) {
    stop("object must be of class 'gcbinary' or 'gctimes'")
  }
  if (!is.null(object$newdata)) {stop("Cannot transport an already transported object")}
  
  fit <- object$calibration$fit 
  model <- object$model
  formula <- object$formula
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
    
    p0 <- c()
    p1 <- c()
    OR <- c()
    delta <- c()
    ratio <- c()
    
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
      
      
      if(model %in% c("lasso","ridge","elasticnet")) {
        .x0 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))),data.valid0)[,-1]
        .x1 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))),data.valid1)[,-1]
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
                outcome=object$outcome,
                group=group,
                n = nrow(newdata),
                nevent = NA,
                p0=p0,
                p1=p1,
                delta=delta,
                ratio=ratio,
                OR=OR,
                effect=effect,
                call=match.call())
    class(res) <- "gcbinary"
    return(res)
  }
  
  
  if (inherits(object, "gctimes")) {
    H0.multi <- object$calibration$H0.multi 
    T.multi <- object$calibration$time
    pro.time <- object$pro.time
    
    AHR <- c()
    RMST0 <- c()
    RMST1 <- c()
    deltaRMST <- c()
    surv0 <- c()
    surv1 <- c()
    deltasurv <- c()
    
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
      

      
      .lp.0 <- tryCatch({
        if (model %in% c("all", "aic", "bic")) { 
          predict(fit, newdata = data.valid0, type="lp")
        } else { 
          .x.valid0 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid0)[,-1] 
          predict(fit, newx = .x.valid0)
        }
      }, error = function(e) {return(NULL)})
      if (is.null(.lp.0)) { BCVerror <- BCVerror + 1 ; next }
      
      .lp.1 <- tryCatch({
        if (model %in% c("all", "aic", "bic")) { 
          predict(fit, newdata = data.valid1, type="lp")
        } else { 
          .x.valid1 <- model.matrix(as.formula(paste("~", deparse(formula[[3]]))), data.valid1)[,-1] 
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
                model=model,
                cv=object$cv, 
                missing=nmiss,
                pro.time=pro.time,
                boot.number = boot.number,
                outcome=object$outcome,
                group=group,
                n = nrow(newdata),
                nevent = NA,
                AHR = AHR,
                RMST0 = RMST0,
                RMST1 = RMST1,
                deltaRMST = deltaRMST,
                surv0 = surv0,
                surv1 = surv1,
                deltasurv = deltasurv,
                call = match.call())
    class(res) <- "gctimes"
    return(res)
  }
}