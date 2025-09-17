gc_survival <- function(formula, data, group, pro.time, effect="ATE", method, param.tune=NULL, cv=10, boot.type="bcv",
                        boot.number=500,  boot.tune=FALSE, progress=TRUE, seed=NULL) {
  # Quality tests
  if(missing(formula)) {stop("The \"formula\" argument is missing (formula)")}
  if(missing(data)) {stop("The \"data\" argument is missing (data.frame)")}
  if(missing(group)) {stop("The \"group\" argument is missing (character string, name of the binary grouping variable in the formula)")}
  if(missing(method)) {stop("Specify one method among : elasticnet, lasso, ridge, all, aic, bic")}
  if(length(method)!=1)
  { stop("Specify one method among : elasticnet, lasso, ridge, all, aic, bic")   }
  
  if(!(method %in% c("elasticnet","lasso","ridge","all","aic","bic"))) {
    stop("Specify one method among : elasticnet, lasso, ridge, all, aic, bic")
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
  formula.all <- formula
  
  datakeep <- data[,which(colnames(data) %in% c(times,failures,group,all_terms))]
  data <- data[,which(colnames(data) %in% c(times,failures,group,all_terms))]
  
  if(is.null(seed)) {seed <- sample(1:1000,1)}
  
  
  if (missing(pro.time) | is.null(pro.time)) {
    warning("The argument \"pro.time\" is missing, the median value of the time variable will be used")
    pro.time <- median(data[,times])
  }
  if (pro.time > max(data[,times])) {stop("The argument \"pro.time\" is higher than the maximum value of the time variable")}
  
  
  
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
  
  if (!(boot.type %in% c("boot","bcv"))) {
    stop("Specify one type to obtain confidence intervals among : boot, bcv")
  }
  
  

  # Initialisation of the tuning parameter


  if(!(method %in% c("lasso","all","ridge","aic","bic","elasticnet"))){
    stop("New \"method\" is not yet implemented, use one among the following : lasso, all, ridge, aic, bic or elasticnet") }
  
  
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
      if (boot.tune & length(param.tune$lambda) == 1) {
        boot.tune <- FALSE
        warning("Only one lambda given, the \"boot.tune\" parameter was set to FALSE")
      }
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
      if (boot.tune & length(param.tune$lambda) == 1) {
        boot.tune <- FALSE
        warning("Only one lambda given, the \"boot.tune\" parameter was set to FALSE")
      }
    }
  }
  
  if (method == "elasticnet") {
    if (!is.list(param.tune) & !is.vector(param.tune) & !is.null(param.tune)) {stop("Tune parameter needs to be a list or a vector (lambda then alpha) or NULL")}
    if (is.list(param.tune) & length(param.tune) != 2) {stop("List tune parameter needs to have a length of 2 (lambda then alpha)")}
    if (is.vector(param.tune) & length(param.tune) != 2) {stop("Vector tune parameter needs to have a length of 2 (lambda then alpha)")}
    if (!is.null(names(param.tune[1]))) {
      if (is.list(param.tune) & names(param.tune[1]) != "lambda" & names(param.tune[2]) != "alpha") {stop("List needs to start with lambda then alpha")}
    }
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
    
    if (boot.tune & length(param.tune$lambda) == 1) {
      boot.tune <- FALSE
      warning("Only one lambda given, the \"boot.tune\" parameter was set to FALSE")
    }
  }

  
  
  ### Effect  

  
  if(effect=="ATE"){ ttt <- which(data[,group] %in% c(0,1))
  }else if(effect=="ATT"){ ttt <- which(data[,group] == 1)
  }else ttt <- which(data[,group] == 0 )
  data <- data[ttt,]
  N <- length(data[,times])
  
  ### Method

  .x <- model.matrix(formula,data)[,-1]
  .y <- Surv(data[,times], data[,failures])
  
  
  set.seed(seed)
  foldid <- sample(rep(seq(cv), length.out = nrow(.x)))
  
if(method == "lasso"){
  if(is.null(param.tune$lambda)==T | length(param.tune$lambda)>1){

    .cv.lasso <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                           nfolds = cv, parallel = FALSE, alpha=1, penalty.factor = .penalty.factor,keep=F,
                           lambda=param.tune$lambda, foldid = foldid)

    .tune.optimal=list(lambda=.cv.lasso$lambda.min)
    rm(.cv.lasso)  }   else{ .tune.optimal=list(lambda=param.tune$lambda) }
}
  if(method == "ridge"){
    if(is.null(param.tune$lambda)==T | length(param.tune$lambda)>1){
      .cv.ridge <- cv.glmnet(x=.x, y=.y, family = "cox",  type.measure = "deviance",
                             parallel = FALSE, alpha=0, penalty.factor = .penalty.factor, nfolds = cv,
                             lambda=param.tune$lambda, foldid = foldid)
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
                                  lambda=param.tune$lambda, foldid = foldid)
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
    formula <- stepAIC( coxph(formula=formula(paste0("Surv(",times,",",failures,")~",group)), data=data),
                     scope=list(lower = formula(paste0("Surv(",times,",",failures,")~",group)), upper = formula),
                     direction="forward", k=2, trace=FALSE)$formula
  } 
  if(method == "bic"){
    formula <- stepAIC( coxph(formula=formula(paste0("Surv(",times,",",failures,")~",group)), data=data),
                        scope=list(lower = formula(paste0("Surv(",times,",",failures,")~",group)), upper = formula),
                        direction="forward", k=log(nrow(data)), trace=FALSE)$formula
  } 

  
  
  ### Calibration survival function
  
  if(method == "all" | method == "aic" | method == "bic") {
    fit <- coxph(formula = formula, data=data)
    
    .lp.coxph <- predict(fit, newdata = data, type="lp")
    .b <- glmnet_basesurv(data[,times], data[,failures], .lp.coxph, centered = FALSE)
    hazard <- .b$cumulative_base_hazard
    fit_times <- .b$times
    .tune.optimal = NULL
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


  if (method == "all" | method == "aic" | method == "bic") {.lp <- predict(fit, newdata = data, type="lp")} else{
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
  
  
  results.surv.calibration <- list(fit=fit, time=T.multi, cumhaz=H.mean, surv=S.mean, H0.multi=H0.multi, lp=lp)
  


  BCVerror <- 0
  pro.time.extrapolate <- 0
  
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
  


    .x.learn = model.matrix(formula.all,data)[,-1][id,]
    .x.valid0 = model.matrix(formula.all,data0)[,-1][-sort(unique(id)),]
    .x.valid1 = model.matrix(formula.all,data1)[,-1][-sort(unique(id)),]
  
  
    
    .y.learn <- Surv(data.learn[,times], data.learn[,failures])

    
    
    ### Unadjusted results
    formula.unadj <- as.formula(paste0("Surv(",times,",",failures,")~",group))
    fit <- coxph(formula = formula.unadj, data=data.learn)
    
    .lp.0.unadj <- predict(fit, newdata = data.valid0, type="lp")
    .lp.1.unadj <- predict(fit, newdata = data.valid1, type="lp")
    
    .lp.unadj_learn <- predict(fit, newdata = data.learn, type="lp")
    .b.unadj <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.unadj_learn, centered = FALSE)
    hazard.unadj <- .b.unadj$cumulative_base_hazard
    fit_times.unadj <- .b.unadj$times
    
    H0.multi.unadj <- c(0, hazard.unadj[fit_times.unadj %in% sort(unique(data.learn[data.learn[,failures]==1,times]))] )
    T.multi.unadj <- c(0, fit_times.unadj[fit_times.unadj %in% sort(unique(data.learn[data.learn[,failures]==1,times]))] )
    
    lp.0.unadj <- as.vector(.lp.0.unadj)
    lp.1.unadj <- as.vector(.lp.1.unadj)
    
    h0.unadj <- (H0.multi.unadj[2:length(T.multi.unadj)] - H0.multi.unadj[1:(length(T.multi.unadj)-1)])
    hi.0.unadj <- exp(lp.0.unadj) * matrix(rep(h0.unadj,length(lp.0.unadj)), nrow=length(lp.0.unadj), byrow=TRUE)
    Si.0.unadj <- exp(-exp(lp.0.unadj) * matrix(rep(H0.multi.unadj,length(lp.0.unadj)), nrow=length(lp.0.unadj), byrow=TRUE))
    hi.0.unadj <- cbind(rep(0,length(lp.0.unadj)),hi.0.unadj)

    h.mean.0.unadj <- apply(Si.0.unadj * hi.0.unadj, FUN="sum", MARGIN=2) / apply(Si.0.unadj, FUN="sum", MARGIN=2)
    H.mean.0.unadj <- cumsum(h.mean.0.unadj)
    S.mean.0.unadj <- exp(-H.mean.0.unadj) 
    
    hi.1.unadj <- exp(lp.1.unadj) * matrix(rep(h0.unadj,length(lp.1.unadj)), nrow=length(lp.1.unadj), byrow=TRUE)
    Si.1.unadj <- exp(-exp(lp.1.unadj) * matrix(rep(H0.multi.unadj,length(lp.1.unadj)), nrow=length(lp.1.unadj), byrow=TRUE))
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
    

    
    ### GC
    if (method == "aic") {
      formula <- stepAIC( coxph(formula=formula(paste0("Surv(",times,",",failures,")~",group)), data=data.learn),
                          scope=list(lower = formula(paste0("Surv(",times,",",failures,")~",group)), upper = formula.all),
                          direction="forward", k=2, trace=FALSE)$formula
      
      fit <- coxph(formula = formula, data=data.learn)
      
      .lp.coxph <- predict(fit, newdata = data.learn, type="lp")
      .b <- glmnet_basesurv(data.learn[,times], data.learn[,failures], .lp.coxph, centered = FALSE)
      hazard <- .b$cumulative_base_hazard
      fit_times <- .b$times
    }
    
    if (method == "bic") {
      formula <- stepAIC( coxph(formula=formula(paste0("Surv(",times,",",failures,")~",group)), data=data.learn),
                          scope=list(lower = formula(paste0("Surv(",times,",",failures,")~",group)), upper = formula.all),
                          direction="forward", k=log(nrow(data.learn)), trace=FALSE)$formula
      
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
        .cv.lasso <- cv.glmnet(x=.x.learn, y=.y.learn, family = "cox",  type.measure = "deviance",
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
        .cv.ridge <- cv.glmnet(x=.x.learn, y=.y.learn, family = "cox",  type.measure = "deviance",
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
          .cv.en<-glmnet::cv.glmnet(x=.x.learn, y=.y.learn, family = "cox",  type.measure = "deviance",
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
    
    if (max(data.learn[,times]) < pro.time) {
      pro.time.extrapolate = pro.time.extrapolate + 1
    }
    
    baseline_hazard <- hazard
    H0.multi <- c(0, baseline_hazard[fit_times %in% sort(unique(data[data[,failures]==1,times]))]  )
    T.multi <- c(0, fit_times[fit_times %in% sort(unique(data[data[,failures]==1,times]))] )
    
  
    
    
    if (method == "all" | method == "aic" | method == "bic") {
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
    
    if (max(T.multi) != max(fit_times)) {
      T.multi = c(T.multi, max(fit_times))
      H.mean.0 = c(H.mean.0, max(H.mean.0))
      S.mean.0 = c(S.mean.0, min(S.mean.0))
      H.mean.1 = c(H.mean.1, max(H.mean.1))
      S.mean.1 = c(S.mean.1, min(S.mean.1))
    }
    
    AHR <- c(AHR, sum(h.mean.1) / sum(h.mean.0) )
    
      .S.mean.0 <- S.mean.0[order(T.multi)]
      .S.mean.1 <- S.mean.1[order(T.multi)]
      .T.multi <- T.multi[order(T.multi)]
      .t <- c(.T.multi[.T.multi <= pro.time], min(pro.time, max(.T.multi)))
      .s0 <- c(.S.mean.0[.T.multi <= pro.time], .S.mean.0[length(.S.mean.0)]) 
      .s1 <- c(.S.mean.1[.T.multi <= pro.time], .S.mean.1[length(.S.mean.1)]) 
      .RMST0 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s0[1:(length(.s0) - 1)])
      .RMST1 <- sum((.t[2:length(.t)] - .t[1:(length(.t) - 1)]) * .s1[1:(length(.s1) - 1)])  
    
      
      
      .surv1 <- .S.mean.1[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
      .surv0 <- .S.mean.0[findInterval(pro.time, T.multi, rightmost.closed = TRUE)]
        
      surv1 <- c(surv1,.surv1)
      surv0 <- c(surv0,.surv0)
      deltasurv <- c(deltasurv, .surv1 - .surv0)
      
    
  RMST0 <- c(RMST0, .RMST0)
  RMST1 <- c(RMST1, .RMST1)

  deltaRMST <- c(deltaRMST, .RMST1 - .RMST0)
    }
    
  if(progress==TRUE){ close(pb) }
  
if (pro.time.extrapolate > 1) {warning(paste0("In at least one boostrap sample the \"pro.time\" was higher than the maximum follow-up time (survival was extrapolated in ",pro.time.extrapolate," bootstrap samples). It is advised to pick a lower value for \"pro.time\""))}
if (BCVerror > 1) {warning(paste0("Skipped ",BCVerror," bootstrap iterations due to the validation dataset containing factors not in the train dataset. Either use type=\"boot\" instead of \"bcv\" or remove factors with rare modalities."))}  
if (!is.null(.warnen)) {warning(paste0("The optimal tuning parameter alpha was equal to ",.warnen,", using ",ifelse(.warnen==0,"ridge","lasso")," instead"))}  
  
  if (method == "aic" | method == "bic") {.tune.optimal = NULL}
  
  

res <- list(calibration=as.list(results.surv.calibration),
            tuning.parameters=.tune.optimal,
            data=datakeep,
            formula=formula.all,
            method=method,
            cv=cv,
            missing=nmiss,
            pro.time=pro.time,
            boot.number = boot.number,
            boot.type = boot.type,
            outcome=list(times=times, failures=failures),
            group=group,
            n = nrow(datakeep),
            nevent = sum(datakeep[,failures]),
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
            b=b,
            call = match.call()
)


class(res) <- "gcsurv"

return(res)
}
