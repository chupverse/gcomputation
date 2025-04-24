
gc_survival <- function(formula, data, group, max.time, effect="ATE", method, param.tune=NULL, cv=10, boot.type=NULL,
                        boot.number=500,  boot.tune=FALSE, progress=TRUE) {
  

  # Quality tests
  
  if(missing(formula)) {stop("The \"formula\" argument is missing (formula)")}
  if(missing(data)) {stop("The \"data\" argument is missing (data.frame)")}
  if(missing(group)) {stop("The \"group\" argument is missing (character string, name of the binary grouping variable in the formula)")}
  if(missing(method)) {stop("Specify one method among : elasticnet, lasso, ridge, all, aic")}
  if(length(method)!=1)
  { stop("Specify one method among : elasticnet, lasso, ridge, all, aic")   }
  
  if(!(method %in% c("elasticnet","lasso","ridge","all","aic"))) {
    stop("Specify one method among : elasticnet, lasso, ridge, all, aic")
  }
  if(!is.data.frame(data)){stop("The argument \"data\" needs to be a data.frame") }
  
  if (!is.logical(boot.tune)) {stop("The argument \"boot.tune\" needs to be a logical TRUE or FALSE")}
  
  mod <- unique(data[,group])
  if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
    stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the \"data\" for argument \"group\"") }
  
  
  if (as.character(class(formula)) != "formula") stop("The argument \"formula\" must be a formula")
  
  times <- as.character(formula[[2]][2]) 
  failures <- as.character(formula[[2]][3])
  all_terms <- attr(terms(formula), "term.labels")
  
  
  if (missing(max.time)) {
    warning("The argument \"max.time\" is missing, the median value of the time variable will be used")
    max.time <- median(data[,times])
  }
  if (max.time > max(data[,times])) {stop("The argument \"max.time\" is higher than the maximum value of the time variable")}
  
  
  
  if(length(grep("tt\\(", all_terms, value = TRUE)) > 0 |
     length(grep("strata\\(", all_terms, value = TRUE)) > 0 |
     length(grep("cluster\\(", all_terms, value = TRUE)) > 0){
    stop("Incorrect argument in the argument \"object\": time-transform
functions, stratification and clustering are not implemented") }
  
  
  
  if( !is.null(group)){
    if(length(group)>1){
      stop("Only one variable can be used as group")
    }
    if(!is.character(group)){
      stop("The argument \"group\" needs to be a character string")
    }
    
    if(!(group %in% colnames(data))){
      stop("Group name is not present in data")
    }
  } else {
    stop("The argument \"group\" needs to be specified")
  }
  if (!(group %in% all_terms)) {
    stop("The argument \"group\" needs to be specified in the formula")
  }
  all_terms <- all_terms[-which(all_terms == group)]
  
  .penalty.factor <- rep(1,length(colnames(model.matrix(formula,data)[,-1])))
  .penalty.factor[which(colnames(model.matrix(formula,data)[,-1]) == group)] <- 0
  
  
  mod <- unique(data[,group])
  if(length(mod) != 2 | ((mod[1] != 0 & mod[2] != 1) & (mod[1] != 1 & mod[2] != 0))){
    stop("Two modalities encoded 0 (for non-treated/non-exposed patients) and 1 (for treated/exposed patients) are required in the \"group\" variable")
  }
  
  
  if(length(data[,times])!=length(data[,failures])){
    stop("The length of the times must be equal to the length of the events in the training data") }
  
  mod2 <- unique(data[,failures])
  if(length(mod2) != 2 | ((mod2[1] != 0 & mod2[2] != 1) & (mod2[1] != 1 & mod2[2] != 0))){
    stop("Two modalities encoded 0 (for censored patients) and 1 (for events) are required in the \"failures\" variable")
  }
  
  if (!is.numeric(data[,times])){
    stop("Time variable is not numeric")}
  
  if (min(data[,times])<=0){
    stop("Time variable needs to be strictly positive")
  }
  
  
  if(length(effect)!=1)
  { stop("Specify one average effect among : ATE, ATT, ATU")   }
  
  if(!(effect %in% c("ATE","ATT","ATU"))) {
    stop("Specify one average effect among : ATE, ATT, ATU")
  }
  
  

  if (any(is.na(data))){
    nmiss <- nrow(data)
    data <- na.omit(data)
    nmiss <- nmiss - nrow(data)
    warning("Rows containing NA values have been removed from the dataset!")
  } else {nmiss <- 0}
  

  if (cv < 3 | !is.numeric(cv)) {
    stop("nfolds must be bigger than 3; nfolds=10 recommended")
  }
  
  if (boot.number < 2 | !is.numeric(boot.number)) {
    stop("boot.number, the number of bootstraps, must be bigger than 2; boot.number=500 recommended")
  }
  
  if (missing(boot.type)) {boot.type = "bcv"}
  if (!(boot.type %in% c("boot","bcv"))) {
    stop("Specify one type to obtain confidence intervals among : boot, bcv")
  }
  
  

  ### Initialisation et recuperation param.tune


  if(!(method %in% c("lasso","all","ridge","aic","elasticnet"))){
    stop("New \"method\" is not yet implemented, use one among the following : lasso, all, ridge, aic or elasticnet") }
  
  
  time.pred <- sort(unique(data[,times]))
  
  if(progress==TRUE){
    max.progess <- boot.number
    pb <- txtProgressBar(min = 0, max = max.progess, style = 3, width = 50, char = "=")
    ip <- 0
    setTxtProgressBar(pb, ip)
  }
  
  
  if(method == "lasso"){
    if(!(is.numeric(param.tune) | is.null(param.tune))){
      stop("Tune parameter lambda for Lasso method needs to be a scalar or a vector or NULL")
    }

    if (is.null(param.tune)) {
      param.tune=list(lambda=NULL)
    } else {
      param.tune = list(lambda = param.tune)
    }
  }  
  
  if(method == "ridge"){
    if(!(is.numeric(param.tune) | is.null(param.tune))){
      stop("Tune parameter lambda for Ridge method needs to be a scalar or a vector or NULL")
    }
    
    if (is.null(param.tune)) {
      param.tune=list(lambda=NULL)
    } else {
      param.tune = list(lambda = param.tune)
    }
  }
  
  if (method == "elasticnet") {
    if (!is.list(param.tune) & !is.vector(param.tune) & !is.null(param.tune)) {stop("Tune parameter needs to be a list or a vector (lambda then alpha) or NULL")}
    if (is.list(param.tune) & length(param.tune) != 2) {stop("List tune parameter needs to have a length of 2 (lambda then alpha)")}
    if (is.vector(param.tune) & length(param.tune) != 2) {stop("Vector tune parameter needs to have a length of 2 (lambda then alpha)")}
    if (is.list(param.tune) & length(param.tune[[1]]) == 1 & length(param.tune[[2]]) > 1) {stop("Lambda needs more than 1 value if more than 1 alpha is provided")} 
    if (is.null(param.tune)) {
      param.tune = list(lambda=NULL, alpha=seq(0,1,.1) )
    } else {
      param.tune = list(lambda=param.tune[[1]], alpha=param.tune[[2]])
    }
    if (!(is.numeric(param.tune$lambda)| is.null(param.tune$lambda))) {
      stop(paste("lambda needs to be a scalar or a vector or NULL"))
    }
    if (!(is.numeric(param.tune$alpha)| is.null(param.tune$alpha))) {
      stop(paste("alpha needs to be a scalar or a vector or NULL"))
    }
  }

  
  
  ### Effect  

  
  if(effect=="ATE"){ ttt <- which(data[,group] %in% c(0,1))
  }else if(effect=="ATT"){ ttt <- which(data[,group] == 1)
  }else ttt <- which(data[,group] == 0 )
  data <- data[ttt,]
  N <- length(data[,times])
  
  ############# Method

  .x <- model.matrix(formula,data)[,-1]
  .y <- Surv(data[,times], data[,failures])
  
if(method == "lasso"){
  if(is.null(param.tune$lambda)==T | length(param.tune$lambda)>1){

    .cv.lasso <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                           nfolds = cv, parallel = FALSE, alpha=1, penalty.factor = .penalty.factor,keep=F,
                           lambda=param.tune$lambda)

    .tune.optimal=list(lambda=.cv.lasso$lambda.min)
    rm(.cv.lasso)  }   else{ .tune.optimal=list(lambda=param.tune$lambda) }
}
  if(method == "ridge"){
    if(is.null(param.tune$lambda)==T | length(param.tune$lambda)>1){
      .cv.ridge <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                             parallel = FALSE, alpha=0, penalty.factor = .penalty.factor, nfolds = cv,
                             lambda=param.tune$lambda)
      .tune.optimal=list(lambda=.cv.ridge$lambda.min)
      rm(.cv.ridge)  } else{ .tune.optimal=list(lambda=param.tune$lambda) }
  }
  .warnen = NULL
  if(method == "elasticnet"){
    if (is.null(param.tune$lambda)==T | length(param.tune$lambda)>1 | length(param.tune$alpha)>1){
      .results<-c()
      for( a in 1:length(param.tune$alpha)){
        .cv.en<-glmnet::cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                                  foldsid="folds", parallel = FALSE, alpha=param.tune$alpha[a],
                                  penalty.factor = .penalty.factor,
                                  lambda=param.tune$lambda)
        .results<-rbind(.results,
                        cbind(rep(param.tune$alpha[a],length(.cv.en$lambda)),.cv.en$lambda,.cv.en$cvm))
      }
      colnames(.results)=c("alpha","lambda","cvm")
      .results=data.frame(.results)
      .tune.optimal=list(alpha=.results[which(.results$cvm==min(.results$cvm)),1][1] ,
                         lambda=.results[which(.results$cvm==min(.results$cvm)),2][1] )
      rm(.cv.en) ; rm(.results) } else{.tune.optimal=list(alpha=param.tune$alpha, lambda=param.tune$lambda) }
    
    
    if (.tune.optimal$alpha == 1 & boot.tune == FALSE) {.warnen=1}
    if (.tune.optimal$alpha == 0 & boot.tune == FALSE) {.warnen=0}
  } 
  if(method == "aic"){
    formula.all = formula
    formula <- stepAIC( coxph(formula=formula(paste0("Surv(",times,",",failures,")~",group)), data=data),
                     scope=list(lower = formula(paste0("Surv(",times,",",failures,")~",group)), upper = formula),
                     direction="forward", k=2, trace=FALSE)$formula
    .tune.optimal = NULL
  } 

  
  
  #### Calibration survival function
  
  if(method == "all" | method == "aic") {
    fit <- coxph(formula = formula, data=data)
    
    .lp.coxph <- predict(fit, newdata = data, type="lp")
    .b <- glmnet_basesurv(data[,times], data[,failures], .lp.coxph, centered = FALSE)
    hazard <- .b$cumulative_base_hazard
    fit_times <- .b$times
  }
  if (method == "lasso") {
    fit <- glmnet(x = .x, y = .y, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                                family = "cox", alpha = 1, penalty.factor = .penalty.factor)
    
    .lp.lasso <- predict(fit, newx = .x)
    .b <- glmnet_basesurv(data[,times], data[,failures], .lp.lasso, centered = FALSE)
    hazard <- .b$cumulative_base_hazard
    fit_times <- .b$times
  }
  if (method == "ridge") {
    fit <- glmnet(x = .x, y = .y, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                  family = "cox", alpha = 0, penalty.factor = .penalty.factor)
    
    .lp.ridge <- predict(fit, newx = .x)
    .b <- glmnet_basesurv(data[,times], data[,failures], .lp.ridge, centered = FALSE)
    hazard <- .b$cumulative_base_hazard
    fit_times <- .b$times
  }
  if (method == "elasticnet") {
    fit <- glmnet(x = .x, y = .y, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                  family = "cox", alpha = .tune.optimal$alpha, penalty.factor = .penalty.factor)
    
    .lp.elasticnet <- predict(fit, newx = .x)
    .b <- glmnet_basesurv(data[,times], data[,failures], .lp.elasticnet, centered = FALSE)
    hazard <- .b$cumulative_base_hazard
    fit_times <- .b$times
  }
  
  
  
  baseline_hazard <- hazard
  H0.multi <- c(0, baseline_hazard[fit_times %in% sort(unique(data[data[,failures]==1,times]))]  )
  T.multi <- c(0, fit_times[fit_times %in% sort(unique(data[data[,failures]==1,times]))] )


  if (method == "all" | method == "aic") {.lp <- predict(fit, newdata = data, type="lp")} else{
    .lp <- predict(fit, newx = .x)
  }
  
  lp <- as.vector(.lp)
  
  h0 <- (H0.multi[2:length(T.multi)] - H0.multi[1:(length(T.multi)-1)]) 
  hi <- exp(lp) * matrix(rep(h0,length(lp)), nrow=length(lp), byrow=TRUE)
  Si <- exp(-exp(lp) * matrix(rep(H0.multi,length(lp)), nrow=length(lp), byrow=TRUE))
  
  hi <- cbind(rep(0,length(lp)),hi)
  
  h.mean <- apply(Si * hi, FUN="sum", MARGIN=2) / apply(Si, FUN="sum", MARGIN=2)
  H.mean <- cumsum(h.mean)
  S.mean <- exp(-H.mean)
  

  if (max(T.multi) != max(fit_times)) {
    T.multi = c(T.multi, max(fit_times))
    H.mean = c(H.mean, max(H.mean))
    S.mean = c(S.mean, min(S.mean))
  }
  
  
  results.surv.calibration <- list(model=fit, time=T.multi, cumhaz=H.mean, surv=S.mean, H0.multi=H0.multi, lp=lp)
  

  ###   Bootstrapping

  delta <- rep(NA, boot.number)
  max.time.extrapolate <- 0
  RMST0 <- c()
  RMST1 <- c()
  deltaRMST <- c() 
  BCVerror <- 0
  for (b in 1:boot.number) {
    
    if(progress == TRUE){
      ip <- ip + 1
      setTxtProgressBar(pb, ip)
    }
    
    id = sample(1:N, size = N, replace = TRUE)
   
    data.learn = data[id,]
    if (boot.type == "bcv") {data.valid = data[-sort(unique(id)),]} else{
      data.valid = data[id,]
    }
    
    if (effect == "ATE")  {
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
    
    data0=data1=data
    data0[,group] = 0
    data1[,group] = 1
  

    ### Fixes the issue when there is modalities not present in valid or not in train
    
    .x.learn = model.matrix(formula,data)[,-1][id,]
    .x.valid0 = model.matrix(formula,data0)[,-1][-sort(unique(id)),]
    .x.valid1 = model.matrix(formula,data1)[,-1][-sort(unique(id)),]
  
  
    
    .y.learn <- Surv(data.learn[,times], data.learn[,failures])

    
    if (method == "aic") {
      formula <- stepAIC( coxph(formula=formula(paste0("Surv(",times,",",failures,")~",group)), data=data),
                          scope=list(lower = formula(paste0("Surv(",times,",",failures,")~",group)), upper = formula.all),
                          direction="forward", k=2, trace=FALSE)$formula
      
      fit <- coxph(formula = formula, data=data.learn)
      
      .lp.coxph <- predict(fit, newdata = data.learn, type="lp")
      .b <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.coxph, centered = FALSE)
      hazard <- .b$cumulative_base_hazard
      fit_times <- .b$times
      
    }
    
    if(method == "all") {
      fit <- coxph(formula = formula, data=data.learn)
      
      .lp.coxph <- predict(fit, newdata = data.learn, type="lp")
      .b <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.coxph, centered = FALSE)
      hazard <- .b$cumulative_base_hazard
      fit_times <- .b$times
    }
    if (method == "lasso") {
      if (boot.tune) {
        .cv.lasso <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                               nfolds = cv, parallel = FALSE, alpha=1, penalty.factor = .penalty.factor,keep=F,
                               lambda=param.tune$lambda)
        .tune.optimal=list(lambda=.cv.lasso$lambda.min)
      }
      fit <- glmnet(x = .x.learn, y = .y.learn, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                    family = "cox", alpha = 1, penalty.factor = .penalty.factor)
      
      
      .lp.lasso <- predict(fit, newx = .x.learn)
      .b <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.lasso, centered = FALSE)
      hazard <- .b$cumulative_base_hazard
      fit_times <- .b$times
    }
    if (method == "ridge") {
      if (boot.tune) {
        .cv.ridge <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                               parallel = FALSE, alpha=0, penalty.factor = .penalty.factor, nfolds = cv,
                               lambda=param.tune$lambda)
        .tune.optimal=list(lambda=.cv.ridge$lambda.min)
      }
      fit <- glmnet(x = .x.learn, y = .y.learn, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                    family = "cox", alpha = 0, penalty.factor = .penalty.factor)
      
      
      .lp.ridge <- predict(fit, newx = .x.learn)
      .b <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.ridge, centered = FALSE)
      hazard <- .b$cumulative_base_hazard
      fit_times <- .b$times
    }
    if (method == "elasticnet") {
      if (boot.tune) {
        .results<-c()
        for( a in 1:length(param.tune$alpha)){
          .cv.en<-glmnet::cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                                    foldsid="folds", parallel = FALSE, alpha=param.tune$alpha[a],
                                    penalty.factor = .penalty.factor,
                                    lambda=param.tune$lambda)
          .results<-rbind(.results,
                          cbind(rep(param.tune$alpha[a],length(.cv.en$lambda)),.cv.en$lambda,.cv.en$cvm))
        }
        colnames(.results)=c("alpha","lambda","cvm")
        .results=data.frame(.results)
        .tune.optimal=list(alpha=.results[which(.results$cvm==min(.results$cvm)),1][1] ,
                           lambda=.results[which(.results$cvm==min(.results$cvm)),2][1] )
      }
      fit <- glmnet(x = .x.learn, y = .y.learn, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                    family = "cox", alpha = .tune.optimal$alpha, penalty.factor = .penalty.factor)
      
      
      .lp.elasticnet <- predict(fit, newx = .x.learn)
      .b <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.elasticnet, centered = FALSE)
      hazard <- .b$cumulative_base_hazard
      fit_times <- .b$times
    }
    
    if (max(data.learn[,times]) < max.time) {
      max.time.extrapolate = max.time.extrapolate + 1
    }
    
    baseline_hazard <- hazard
    H0.multi <- c(0, baseline_hazard[fit_times %in% sort(unique(data[data[,failures]==1,times]))]  )
    T.multi <- c(0, fit_times[fit_times %in% sort(unique(data[data[,failures]==1,times]))] )
    
  
    
    
    if (method == "all" | method == "aic") {
      .lp.0 <- tryCatch({predict(fit, newdata = data.valid0, type="lp")}, error = function(e) {return(NULL) })
      if (is.null(.lp.0)) {BCVerror <- BCVerror + 1 ; next}
      .lp.1 <- predict(fit, newdata = data.valid1, type="lp")
    } else{ 
      .lp.0 <- predict(fit, newx = .x.valid0) 
      .lp.1 <- predict(fit, newx = .x.valid1)
    }
    
    lp.0 <- as.vector(.lp.0)
    lp.1 <- as.vector(.lp.1)
    
    h0 <- (H0.multi[2:length(T.multi)] - H0.multi[1:(length(T.multi)-1)]) 
    
    hi.0 <- exp(lp.0) * matrix(rep(h0,length(lp.0)), nrow=length(lp.0), byrow=TRUE)
    Si.0 <- exp(-exp(lp.0) * matrix(rep(H0.multi,length(lp.0)), nrow=length(lp.0), byrow=TRUE))
    hi.0 <- cbind(rep(0,length(lp.0)),hi.0)
    hi.1 <- exp(lp.1) * matrix(rep(h0,length(lp.1)), nrow=length(lp.1), byrow=TRUE)
    Si.1 <- exp(-exp(lp.1) * matrix(rep(H0.multi,length(lp.1)), nrow=length(lp.1), byrow=TRUE))
    hi.1 <- cbind(rep(0,length(lp.1)),hi.1)
    
    h.mean.0 <- apply(Si.0 * hi.0, FUN="sum", MARGIN=2) / apply(Si.0, FUN="sum", MARGIN=2)
    H.mean.0 <- cumsum(h.mean.0)
    S.mean.0 <- exp(-H.mean.0)
    h.mean.1 <- apply(Si.1 * hi.1, FUN="sum", MARGIN=2) / apply(Si.1, FUN="sum", MARGIN=2)
    H.mean.1 <- cumsum(h.mean.1)
    S.mean.1 <- exp(-H.mean.1)
    
    ### Extrapolate max time
    if (max(T.multi) != max(fit_times)) {
      T.multi = c(T.multi, max(fit_times))
      H.mean.0 = c(H.mean.0, max(H.mean.0))
      S.mean.0 = c(S.mean.0, min(S.mean.0))
      H.mean.1 = c(H.mean.1, max(H.mean.1))
      S.mean.1 = c(S.mean.1, min(S.mean.1))
    }
    

      .S.mean.0 <- S.mean.0[order(T.multi)]
      .S.mean.1 <- S.mean.1[order(T.multi)]
      .T.multi <- T.multi[order(T.multi)]
      .t <- c(.T.multi[.T.multi <= max.time], min(max.time, max(.T.multi)))
      .s0 <- c(.S.mean.0[.T.multi <= max.time], .S.mean.0[length(.S.mean.0)]) 
      .s1 <- c(.S.mean.1[.T.multi <= max.time], .S.mean.1[length(.S.mean.1)]) 
      .RMST0 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s0[1:(length(.s0) - 1)])
      .RMST1 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s1[1:(length(.s1) - 1)])  
    
    
  RMST0 <- c(RMST0, .RMST0)
  RMST1 <- c(RMST1, .RMST1)

  deltaRMST <- c(deltaRMST, .RMST1 - .RMST0)
    }
    
  if(progress==TRUE){ close(pb) }

if (max.time.extrapolate > 1) {warning(paste0("In at least one boostrap sample the \"max.time\" was higher than the maximum follow-up time (survival was extrapolated in ",max.time.extrapolate," bootstrap samples). It is advised to pick a lower value for \"max.time\""))}
if (BCVerror > 1) {warning(paste0("Skipped ",BCVerror," bootstrap iterations due to the validation dataset containing factors not in the train dataset. Either use type=\"boot\" instead of \"bcv\" or remove factors with rare modalities."))}  
if (!is.null(.warnen)) {warning(paste0("The optimal tuning parameter alpha was equal to ",.warnen,", using ",ifelse(.warnen==0,"ridge","lasso")," instead"))}  
  
  if (method == "aic") {
    .tune.optimal = NULL
    formula = formula.all
  }
  
res <- list(calibration=as.list(as.list(results.surv.calibration)),
            tuning.parameters=.tune.optimal,
            data=data[,which(colnames(data) %in% c(times,failures,group,all_terms))],
            formula=formula,
            method=method,
            cv=cv,
            missing=nmiss,
            max.time=max.time,
            boot.number = boot.number,
            outcomes=list(times=times, failures=failures),
            RMST = list(all.RMST0 = RMST0,
                        mean.RMST0=mean(RMST0, na.rm=TRUE),
                        se.RMST0 = sd(RMST0, na.rm=TRUE),
                        ci.low.asympt.RMST0 = mean(RMST0, na.rm=TRUE) - qnorm(0.975, 0, 1)*sd(RMST0, na.rm=TRUE),
                        ci.upp.asympt.RMST0 = mean(RMST0, na.rm=TRUE) + qnorm(0.975, 0, 1)*sd(RMST0, na.rm=TRUE),
                        ci.low.nonpara.RMST0 = quantile(RMST0, probs = 0.025, na.rm = T),
                        ci.upp.nonpara.RMST0 = quantile(RMST0, probs = 0.975, na.rm = T),
                        
                        all.RMST1 = RMST1,
                        mean.RMST1=mean(RMST1, na.rm=TRUE),
                        se.RMST1 = sd(RMST1, na.rm=TRUE),
                        ci.low.asympt.RMST1 = mean(RMST1, na.rm=TRUE) - qnorm(0.975, 0, 1)*sd(RMST1, na.rm=TRUE),
                        ci.upp.asympt.RMST1 = mean(RMST1, na.rm=TRUE) + qnorm(0.975, 0, 1)*sd(RMST1, na.rm=TRUE),
                        ci.low.nonpara.RMST1 = quantile(RMST1, probs = 0.025, na.rm = T),
                        ci.upp.nonpara.RMST1 = quantile(RMST1, probs = 0.975, na.rm = T),
              
                        all.deltaRMST=deltaRMST, 
                        mean.deltaRMST=mean(deltaRMST, na.rm=TRUE),
                         se.deltaRMST = sd(deltaRMST, na.rm=TRUE),
                         ci.low.asympt.deltaRMST = mean(deltaRMST, na.rm=TRUE) - qnorm(0.975, 0, 1)*sd(deltaRMST, na.rm=TRUE),
                         ci.upp.asympt.deltaRMST = mean(deltaRMST, na.rm=TRUE) + qnorm(0.975, 0, 1)*sd(deltaRMST, na.rm=TRUE),
                        ci.low.nonpara.deltaRMST = quantile(deltaRMST, probs = 0.025, na.rm = T),
                        ci.upp.nonpara.deltaRMST = quantile(deltaRMST, probs = 0.975, na.rm = T),
                         p.value.deltaRMST = ifelse(mean(deltaRMST, na.rm=TRUE)/sd(deltaRMST, na.rm=TRUE)<0,2*pnorm(mean(deltaRMST, na.rm=TRUE)/sd(deltaRMST, na.rm=TRUE)),2*(1-pnorm(mean(deltaRMST, na.rm=TRUE)/sd(deltaRMST, na.rm=TRUE))))
                          )
            )

class(res) <- "gcsurv"

return(res)
}

