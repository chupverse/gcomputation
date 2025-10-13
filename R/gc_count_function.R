#' Title
#'
#' @param formula 
#' @param data 
#' @param group 
#' @param effect 
#' @param model 
#' @param param.tune 
#' @param cv 
#' @param boot.type 
#' @param boot.number 
#' @param boot.tune 
#' @param progress 
#' @param seed 
#'
#' @returns
#' @export
#'
#' @examples
.gc_count <- function(formula, data, group, effect="ATE", model, param.tune=NULL, cv=10, boot.type="bcv",
                           boot.number=500, boot.tune=FALSE, progress=TRUE, seed=NULL) {
  # Quality tests
  if(missing(formula)) {stop("The \"formula\" argument is missing (formula)")}
  if(missing(data)) {stop("The \"data\" argument is missing (data.frame)")}
  if(missing(group)) {stop("The \"group\" argument is missing (character string, name of the binary grouping variable in the formula)")}
  if(missing(model)) {stop("Specify one model among : elasticnet, lasso, ridge, all, aic, bic")}
  if(length(model)!=1)
  { stop("Specify one model among : elasticnet, lasso, ridge, all, aic, bic")   }
  
  if(!(model %in% c("elasticnet","lasso","ridge","all","aic", "bic"))) {
    stop("Specify one model among : elasticnet, lasso, ridge, all, aic, bic")
  }
  if(!is.data.frame(data)){stop("The argument \"data\" needs to be a data.frame") }
  
  if (!is.logical(boot.tune)) {stop("The argument \"boot.tune\" needs to be a logical TRUE or FALSE")}
  
  
  if (as.character(class(formula)) != "formula") stop("The argument \"formula\" must be a formula")
  
  outcome <- as.character(formula[[2]])
  all_terms <- attr(terms(formula), "term.labels")
  formula.all <- formula
  
  datakeep <- data[,which(colnames(data) %in% all.vars(formula))]
  data <- data[,which(colnames(data) %in% all.vars(formula))]
  
  if(is.null(seed)) {seed <- sample(1:1000,1)}
  
  
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
  
  # verification that the outcome is numeric
  
  if(!is.numeric(data[,outcome])){
    stop("The outcome variable needs to be numeric")
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
  
  
  if(!(model %in% c("lasso","all","ridge","aic","bic","elasticnet"))){
    stop("New \"model\" is not yet implemented, use one among the following : lasso, all, ridge, aic, bic or elasticnet") }
  
  
  
  if(progress==TRUE){
    max.progess <- boot.number
    pb <- txtProgressBar(min = 0, max = max.progess, style = 3, width = 50, char = "=")
    ip <- 0
    setTxtProgressBar(pb, ip)
  }
  
  
  
  
  if(model == "lasso"){
    if(!(is.numeric(param.tune) | is.null(param.tune))){
      stop("Tune parameter lambda for Lasso model needs to be a scalar or a vector or NULL")
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
  
  if(model == "ridge"){
    if(!(is.numeric(param.tune) | is.null(param.tune))){
      stop("Tune parameter lambda for Ridge model needs to be a scalar or a vector or NULL")
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
  
  if (model == "elasticnet") {
    if (!is.list(param.tune) & !is.vector(param.tune) & !is.null(param.tune)) {stop("Tune parameter needs to be a list or a vector (lambda then alpha) or NULL")}
    if (is.list(param.tune) & length(param.tune) != 2) {stop("List tune parameter needs to have a length of 2 (lambda then alpha)")}
    if (is.vector(param.tune) & length(param.tune) != 2) {stop("Vector tune parameter needs to have a length of 2 (lambda then alpha)")}
    if (!is.null(names(param.tune[1]))) {
      if (is.list(param.tune) & (names(param.tune[1]) != "lambda" | names(param.tune[2]) != "alpha")) {stop("List needs to start with lambda then alpha")}
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
  N <- length(data[,outcome])
  
  ### model
  
  .x <- model.matrix(formula,data)[,-1]
  .y <- data[,outcome]
  
  
  set.seed(seed)
  foldid <- sample(rep(seq(cv), length.out = nrow(.x)))
  
  if(model == "lasso"){
    if(is.null(param.tune$lambda)==T | length(param.tune$lambda)>1){
      
      .cv.lasso <- cv.glmnet(x=.x, y=.y, family = "poisson",  type.measure = "deviance",
                             nfolds = cv, parallel = FALSE, alpha=1, penalty.factor = .penalty.factor,keep=F,
                             lambda=param.tune$lambda, foldid = foldid)
      
      .tune.optimal=list(lambda=.cv.lasso$lambda.min)
      rm(.cv.lasso)  }   else{ .tune.optimal=list(lambda=param.tune$lambda) }
  }
  if(model == "ridge"){
    if(is.null(param.tune$lambda)==T | length(param.tune$lambda)>1){
      .cv.ridge <- cv.glmnet(x=.x, y=.y, family = "poisson",  type.measure = "deviance",
                             parallel = FALSE, alpha=0, penalty.factor = .penalty.factor, nfolds = cv,
                             lambda=param.tune$lambda, foldid = foldid)
      .tune.optimal=list(lambda=.cv.ridge$lambda.min)
      rm(.cv.ridge)  } else{ .tune.optimal=list(lambda=param.tune$lambda) }
  }
  .warnen = NULL
  if(model == "elasticnet"){
    if (is.null(param.tune$lambda)==T | length(param.tune$lambda)>1 | length(param.tune$alpha)>1){
      .results<-c()
      for( a in 1:length(param.tune$alpha)){
        .cv.en<-glmnet::cv.glmnet(x=.x, y=.y, family = "poisson",  type.measure = "deviance",
                                  parallel = FALSE, alpha=param.tune$alpha[a],
                                  penalty.factor = .penalty.factor,
                                  lambda=param.tune$lambda, foldid = foldid)
        .results<-rbind(.results,
                        cbind(rep(param.tune$alpha[a],length(.cv.en$lambda)),.cv.en$lambda,.cv.en$cvm))
      }
      colnames(.results)=c("alpha","lambda","cvm")
      .results=data.frame(.results)
      .tune.optimal=list(alpha=.results[which(.results$cvm==min(.results$cvm)),1][1] ,
                         lambda=.results[which(.results$cvm==min(.results$cvm)),2][1] )
      rm(.cv.en) ; rm(.results) } else{.tune.optimal=list(alpha=param.tune$alpha, lambda=param.tune$lambda)}
    
    
    if (.tune.optimal$alpha == 1 & boot.tune == FALSE) {.warnen=1}
    if (.tune.optimal$alpha == 0 & boot.tune == FALSE) {.warnen=0}
  } 
  
  if(model == "aic"){
    formula <- stepAIC( glm(formula=formula(paste0(outcome,"~",group)), data=data,family="poisson"),
                        scope=list(lower = formula(paste0(outcome,"~",group)), upper = formula),
                        direction="forward", k=2, trace=FALSE)$formula
  } 
  if(model == "bic"){
    formula <- stepAIC( glm(formula=formula(paste0(outcome,"~",group)), data=data,family="poisson"),
                        scope=list(lower = formula(paste0(outcome,"~",group)), upper = formula),
                        direction="forward", k=log(nrow(data)), trace=FALSE)$formula
  } 
  
  
  
  #### Calibration linear function
  
  
  if(model == "all" | model == "aic" | model == "bic") {
    .tune.optimal = formula
    fit <- glm(formula = formula, data=data, family="poisson")
    calibration.predict <- predict(fit, newdata = data, type = "response")
  }
  if (model == "lasso") {
    fit <- glmnet(x = .x, y = .y, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                  family = "poisson", alpha = 1, penalty.factor = .penalty.factor)
    calibration.predict <- predict(fit, newx=.x, type="response")
  }
  if (model == "ridge") {
    fit <- glmnet(x = .x, y = .y, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                  family = "poisson", alpha = 0, penalty.factor = .penalty.factor)
    calibration.predict <- predict(fit, newx=.x, type="response")
  }
  if (model == "elasticnet") {
    fit <- glmnet(x = .x, y = .y, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                  family = "poisson", alpha = .tune.optimal$alpha, penalty.factor = .penalty.factor)
    calibration.predict <- predict(fit, newx=.x, type="response")
  }
  
  .tune.optimal.totalpop <- .tune.optimal
  calibration.fit <- fit
  
  
  ###   Bootstrapping
  
  BCVerror <- 0
  c0 <- c()
  c1 <- c()
  
  delta <- c()
  ratio <- c()
  
  c0.unadj <- c()
  c1.unadj <- c()
  
  delta.unadj <- c()
  ratio.unadj <- c()
  
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
    
    if (boot.type == "bcv") {
      .x.learn  = model.matrix(formula.all, data)[,-1][id,]
      .x.valid0 = model.matrix(formula.all, data0)[,-1][-sort(unique(id)),]
      .x.valid1 = model.matrix(formula.all, data1)[,-1][-sort(unique(id)),]
      
      if (ncol(model.matrix(formula.all, droplevels(data.learn))) < ncol(model.matrix(formula.all, data.valid))) {BCVerror <- BCVerror + 1 ; next} 
    } else {
      .x.learn  = model.matrix(formula.all, data)[,-1][id,]
      .x.valid0 = model.matrix(formula.all, data0)[,-1][id,]
      .x.valid1 = model.matrix(formula.all, data1)[,-1][id,]
    }
    
    
    .y.learn <- data.learn[,outcome]
    
    
    
    
    ### Unadjusted results
    fit <- glm(formula = as.formula(paste(outcome,"~",group)), data=data.learn, family="poisson")
    
    .c0.unadj = mean(predict(fit, newdata = data.valid0, type = "response"))
    .c1.unadj = mean(predict(fit, newdata = data.valid1, type = "response"))
    
    .delta.unadj = .c1.unadj - .c0.unadj
    .ratio.unadj = .c1.unadj / .c0.unadj
    
    
    ### GC
    if (model == "aic") {
      formula <- stepAIC( glm(formula=formula(paste0(outcome,"~",group)), data=data.learn,family="poisson"),
                          scope=list(lower = formula(paste0(outcome,"~",group)), upper = formula.all),
                          direction="forward", k=2, trace=FALSE)$formula
      
      fit <- suppressWarnings(glm(formula = formula, data=data.learn, family="poisson"))
      
      .c0 = suppressWarnings(mean(predict(fit, newdata = data.valid0, type = "response")))
      .c1 = suppressWarnings(mean(predict(fit, newdata = data.valid1, type = "response")))
      
      .delta = .c1 - .c0
      .ratio = .c1 / .c0
    }
    
    if (model == "bic") {
      formula <- stepAIC( glm(formula=formula(paste0(outcome,"~",group)), data=data.learn,family="poisson"),
                          scope=list(lower = formula(paste0(outcome,"~",group)), upper = formula.all),
                          direction="forward", k=log(nrow(data.learn)), trace=FALSE)$formula
      
      fit <- suppressWarnings(glm(formula = formula, data=data.learn, family="poisson"))
      
      .c0 = suppressWarnings(mean(predict(fit, newdata = data.valid0, type = "response")))
      .c1 = suppressWarnings(mean(predict(fit, newdata = data.valid1, type = "response")))
      
      .delta = .c1 - .c0
      .ratio = .c1 / .c0
    }
    
    if(model == "all") {
      fit <- suppressWarnings(glm(formula = formula, data=data.learn, family="poisson"))
      
      .c0 = suppressWarnings(mean(predict(fit, newdata = data.valid0, type = "response")))
      .c1 = suppressWarnings(mean(predict(fit, newdata = data.valid1, type = "response")))
      
      .delta = .c1 - .c0
      .ratio = .c1 / .c0
    }
    
    if (model == "lasso") {
      if (boot.tune) {
        .cv.lasso <- cv.glmnet(x=.x.learn, y=.y.learn, family = "poisson",  type.measure = "deviance",
                               nfolds = cv, parallel = FALSE, alpha=1, penalty.factor = .penalty.factor,keep=F,
                               lambda=param.tune$lambda)
        .tune.optimal=list(lambda=.cv.lasso$lambda.min)
      }
      fit <- glmnet(x = .x.learn, y = .y.learn, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                    family = "poisson", alpha = 1, penalty.factor = .penalty.factor)
      
      .c0 = mean(predict(fit, newx=.x.valid0, type="response"))
      .c1 = mean(predict(fit, newx=.x.valid1, type="response"))
      
      .delta = .c1 - .c0
      .ratio = .c1 / .c0
    }
    
    if (model == "ridge") {
      if (boot.tune) {
        .cv.ridge <- cv.glmnet(x=.x.learn, y=.y.learn, family = "poisson",  type.measure = "deviance",
                               parallel = FALSE, alpha=0, penalty.factor = .penalty.factor, nfolds = cv,
                               lambda=param.tune$lambda)
        .tune.optimal=list(lambda=.cv.ridge$lambda.min)
      }
      fit <- glmnet(x = .x.learn, y = .y.learn, lambda = .tune.optimal$lambda,  type.measure = "deviance",
                    family = "poisson", alpha = 0, penalty.factor = .penalty.factor)
      
      .c0 = mean(predict(fit, newx=.x.valid0, type="response"))
      .c1 = mean(predict(fit, newx=.x.valid1, type="response"))
      
      .delta = .c1 - .c0
      .ratio = .c1 / .c0
    }
    if (model == "elasticnet") {
      if (boot.tune) {
        .results<-c()
        for( a in 1:length(param.tune$alpha)){
          .cv.en<-glmnet::cv.glmnet(x=.x.learn, y=.y.learn, family = "poisson",  type.measure = "deviance",
                                    parallel = FALSE, alpha=param.tune$alpha[a],
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
                    family = "poisson", alpha = .tune.optimal$alpha, penalty.factor = .penalty.factor)
      
      .c0 = mean(predict(fit, newx=.x.valid0, type="response"))
      .c1 = mean(predict(fit, newx=.x.valid1, type="response"))
      
      .delta = .c1 - .c0
      .ratio = .c1 / .c0
    }
    
    c0 <- c(c0, .c0)
    c1 <- c(c1, .c1)
    
    delta <- c(delta, .delta)
    ratio <- c(ratio, .ratio)
    
    c0.unadj <- c(c0.unadj, .c0.unadj)
    c1.unadj <- c(c1.unadj, .c1.unadj)
    
    delta.unadj <- c(delta.unadj, .delta.unadj)
    ratio.unadj <- c(ratio.unadj, .ratio.unadj)
  }
  
  if(progress==TRUE){ close(pb) }
  
  if (BCVerror > 1) {warning(paste0("Skipped ",BCVerror," bootstrap iterations and only used ", boot.number-BCVerror," iterations due to the validation dataset containing factors not in the train dataset. Either use type=\"boot\" instead of \"bcv\" or remove factors with rare modalities."))}  
  if (!is.null(.warnen)) {warning(paste0("The optimal tuning parameter alpha was equal to ",.warnen,", using ",ifelse(.warnen==0,"ridge","lasso")," instead"))}  
  
  
  
  res <- list(calibration=list(fit=calibration.fit,predict=calibration.predict),
              tuning.parameters=.tune.optimal.totalpop,
              data=datakeep,
              formula=formula.all,
              model=model,
              cv=cv,
              missing=nmiss,
              boot.number = boot.number,
              boot.type = boot.type,
              outcome=outcome,
              group=group,
              n = nrow(datakeep),
              mean_outcome = mean(datakeep[,outcome]),
              c0 = c0,  
              c1 = c1,
              delta = delta,
              ratio = ratio,
              c0.unadj = c0.unadj,
              c1.unadj = c1.unadj,
              delta.unadj = delta.unadj,
              ratio.unadj = ratio.unadj,
              call = match.call()
  )
  
  class(res) <- "gccount"
  
  return(res)
}
