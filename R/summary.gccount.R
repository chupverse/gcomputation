summary.gccount <- function (object, digits=4, ci.type=NULL, ci.level=0.95, unadjusted=TRUE, ...)
{
  if (!is.null(ci.type)){
    if (length(ci.type) != 1) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
    if (!(ci.type %in% c("norm","perc"))) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
  }
  
  x <- object
  if (x$model %in% c("lasso","ridge","elasticnet")) {cat("model: ",x$model,", tuning parameters: ", sep = "")
    if (x$model == "elasticnet") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep="")}
    if (x$model == "lasso") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep="")}
    if (x$model == "ridge") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep="")}
    cat("\nCall:", "\n", sep = "")
    dput(x$formula)}
  if (x$model %in% c("all","aic","bic")) {cat(x$model," model \nCall:", "\n", sep = "")
    dput(x$tuning.parameters)}
  cat("\n")
  
  
  cat("G-computation : \n")
  tmp <- matrix(c(mean(x$c0, na.rm=TRUE), sd(x$c0, na.rm=TRUE), mean(x$c0, na.rm=TRUE)/sd(x$c0, na.rm=TRUE), NA), nrow=1)
  colnames(tmp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp) <- "c0"
  res <- tmp
  
  tmp <- matrix(c(mean(x$c1, na.rm=TRUE), sd(x$c1, na.rm=TRUE), mean(x$c1, na.rm=TRUE)/sd(x$c1, na.rm=TRUE), NA), nrow=1)
  rownames(tmp) <- "c1"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$delta, na.rm=TRUE), sd(x$delta, na.rm=TRUE), mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE),
                  ifelse(mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE)<0,
                         2*pnorm(mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "c1-c0"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$ratio, na.rm=TRUE), sd(x$ratio, na.rm=TRUE), mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE),
                  ifelse(mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE)<0,
                         2*pnorm(mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "c1/c0"
  res <- rbind(res,tmp)
  
  
  if (!is.null(ci.type)) {
    if (ci.type == "norm") {
      ci_vals <- matrix(c(mean(x$c0, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c0, na.rm=TRUE),
                          mean(x$c0, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c0, na.rm=TRUE),
                          mean(x$c1, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c1, na.rm=TRUE),
                          mean(x$c1, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c1, na.rm=TRUE),
                          mean(x$delta, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta, na.rm=TRUE),
                          mean(x$delta, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta, na.rm=TRUE),
                          mean(x$ratio, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio, na.rm=TRUE),
                          mean(x$ratio, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio, na.rm=TRUE)
      ), ncol=2, byrow=TRUE)
    }
    if (ci.type == "perc") {
      ci_vals <- matrix(c(quantile(x$c0, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$c0, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$c1, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$c1, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$delta, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$delta, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$ratio, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$ratio, probs = 1-(1-ci.level)/2, na.rm = TRUE)
      ), ncol=2, byrow=TRUE)
    }
    tmp <- cbind(res[,1], `Lower CI` = ci_vals[,1], `Upper CI` = ci_vals[,2])
    colnames(tmp)[1] <- "Estimate"
    res <- cbind(tmp, res[,2:4])
  }
  
  printCoefmat(res, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  cat("\n")
  
  res_GC <- res
  
  if (!is.null(object$newdata)) {unadjusted = FALSE}
  if (unadjusted == TRUE) {
    cat("Unadjusted : \n")
    tmp <- matrix(c(mean(x$c0.unadj, na.rm=TRUE), sd(x$c0.unadj, na.rm=TRUE), mean(x$c0.unadj, na.rm=TRUE)/sd(x$c0.unadj, na.rm=TRUE), NA), nrow=1)
    colnames(tmp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp) <- "c0"
    res <- tmp
    
    tmp <- matrix(c(mean(x$c1.unadj, na.rm=TRUE), sd(x$c1.unadj, na.rm=TRUE), mean(x$c1.unadj, na.rm=TRUE)/sd(x$c1.unadj, na.rm=TRUE), NA), nrow=1)
    rownames(tmp) <- "c1"
    res <- rbind(res,tmp)
    
    tmp <- matrix(c(mean(x$delta.unadj, na.rm=TRUE), sd(x$delta.unadj, na.rm=TRUE), mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE),
                    ifelse(mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE)<0,
                           2*pnorm(mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE)),
                           2*(1-pnorm(mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE))))), nrow=1)
    rownames(tmp) <- "c1-c0"
    res <- rbind(res,tmp)
    
    tmp <- matrix(c(mean(x$ratio.unadj, na.rm=TRUE), sd(x$ratio.unadj, na.rm=TRUE), mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE),
                    ifelse(mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE)<0,
                           2*pnorm(mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE)),
                           2*(1-pnorm(mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE))))), nrow=1)
    rownames(tmp) <- "c1/c0"
    res <- rbind(res,tmp)
    
    if (!is.null(ci.type)) {
      if (ci.type == "norm") {
        ci_vals <- matrix(c(mean(x$c0.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c0.unadj, na.rm=TRUE),
                            mean(x$c0.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c0.unadj, na.rm=TRUE),
                            mean(x$c1.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c1.unadj, na.rm=TRUE),
                            mean(x$c1.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$c1.unadj, na.rm=TRUE),
                            mean(x$delta.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta.unadj, na.rm=TRUE),
                            mean(x$delta.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta.unadj, na.rm=TRUE),
                            mean(x$ratio.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio.unadj, na.rm=TRUE),
                            mean(x$ratio.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio.unadj, na.rm=TRUE)
        ), ncol=2, byrow=TRUE)
      }
      if (ci.type == "perc") {
        ci_vals <- matrix(c(quantile(x$c0.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$c0.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                            quantile(x$c1.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$c1.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                            quantile(x$delta.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$delta.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                            quantile(x$ratio.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$ratio.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE)
        ), ncol=2, byrow=TRUE)
      }
      tmp <- cbind(res[,1], `Lower CI` = ci_vals[,1], `Upper CI` = ci_vals[,2])
      colnames(tmp)[1] <- "Estimate"
      res <- cbind(tmp, res[,2:4])
    }
    
    printCoefmat(res, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  }
  
  
  cat("\n")
  if (!is.na(object$mean_outcome)) {
    cat(paste0("n= ",x$n,", mean of the outcome= ",x$mean_outcome))
  } else {
    cat(paste0("n= ",x$n))
  }
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep=""); cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep=""); cat("\n") }
  
  if (unadjusted == TRUE) {
    invisible(list(
      GC = res_GC,
      unadjusted = res))
  } else {
    invisible(list(
      GC = res_GC))
  }
}
