summary.gcbinary <- function (object, digits=4, ci.type=NULL, ci.level=0.95, unadjusted=TRUE, ...)
{
  if (!is.null(ci.type)){
    if (length(ci.type) != 1) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
    if (!(ci.type %in% c("norm","perc"))) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
  }

  x <- object
  if (x$model %in% c("lasso","ridge","elasticnet")) {cat("model: ",x$model,", tuning parameters: ", sep = "")
    if (x$model == "elasticnet") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep=" ")}
    if (x$model == "lasso") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep=" ")}
    if (x$model == "ridge") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep=" ")}
    cat("\nCall:", "\n", sep = "")
    dput(x$formula)}
  if (x$model %in% c("all","aic","bic")) {cat(x$model," model \nCall:", "\n", sep = "")
    dput(x$tuning.parameters)}
  cat("\n")

  
  cat("G-computation : \n")
  tmp <- matrix(c(mean(x$adjusted.results$p0, na.rm=TRUE), sd(x$adjusted.results$p0, na.rm=TRUE), mean(x$adjusted.results$p0, na.rm=TRUE)/sd(x$adjusted.results$p0, na.rm=TRUE), NA), nrow=1)
  colnames(tmp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp) <- "P0"
  res <- tmp
  
  tmp <- matrix(c(mean(x$adjusted.results$p1, na.rm=TRUE), sd(x$adjusted.results$p1, na.rm=TRUE), mean(x$adjusted.results$p1, na.rm=TRUE)/sd(x$adjusted.results$p1, na.rm=TRUE), NA), nrow=1)
  rownames(tmp) <- "P1"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$adjusted.results$delta, na.rm=TRUE), sd(x$adjusted.results$delta, na.rm=TRUE), mean(x$adjusted.results$delta, na.rm=TRUE)/sd(x$adjusted.results$delta, na.rm=TRUE),
                  ifelse(mean(x$adjusted.results$delta, na.rm=TRUE)/sd(x$adjusted.results$delta, na.rm=TRUE)<0,
                         2*pnorm(mean(x$adjusted.results$delta, na.rm=TRUE)/sd(x$adjusted.results$delta, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$adjusted.results$delta, na.rm=TRUE)/sd(x$adjusted.results$delta, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "P1-P0"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$adjusted.results$ratio, na.rm=TRUE), sd(x$adjusted.results$ratio, na.rm=TRUE), mean(x$adjusted.results$ratio, na.rm=TRUE)/sd(x$adjusted.results$ratio, na.rm=TRUE),
                  ifelse(mean(x$adjusted.results$ratio, na.rm=TRUE)/sd(x$adjusted.results$ratio, na.rm=TRUE)<0,
                         2*pnorm(mean(x$adjusted.results$ratio, na.rm=TRUE)/sd(x$adjusted.results$ratio, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$adjusted.results$ratio, na.rm=TRUE)/sd(x$adjusted.results$ratio, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "P1/P0"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$adjusted.results$OR, na.rm=TRUE), sd(x$adjusted.results$OR, na.rm=TRUE), mean(x$adjusted.results$OR, na.rm=TRUE)/sd(x$adjusted.results$OR, na.rm=TRUE),
                  ifelse(mean(x$adjusted.results$OR, na.rm=TRUE)/sd(x$adjusted.results$OR, na.rm=TRUE)<0,
                         2*pnorm(mean(x$adjusted.results$OR, na.rm=TRUE)/sd(x$adjusted.results$OR, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$adjusted.results$OR, na.rm=TRUE)/sd(x$adjusted.results$OR, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "OR"
  res <- rbind(res,tmp)
  
  
  if (!is.null(ci.type)) {
    if (ci.type == "norm") {
      ci_vals <- matrix(c(mean(x$adjusted.results$p0, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$p0, na.rm=TRUE),
                          mean(x$adjusted.results$p0, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$p0, na.rm=TRUE),
                          mean(x$adjusted.results$p1, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$p1, na.rm=TRUE),
                          mean(x$adjusted.results$p1, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$p1, na.rm=TRUE),
                          mean(x$adjusted.results$delta, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$delta, na.rm=TRUE),
                          mean(x$adjusted.results$delta, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$delta, na.rm=TRUE),
                          mean(x$adjusted.results$ratio, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$ratio, na.rm=TRUE),
                          mean(x$adjusted.results$ratio, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$ratio, na.rm=TRUE),
                          mean(x$adjusted.results$OR, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$OR, na.rm=TRUE),
                          mean(x$adjusted.results$OR, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$adjusted.results$OR, na.rm=TRUE)
                          ), ncol=2, byrow=TRUE)
    }
    if (ci.type == "perc") {
      ci_vals <- matrix(c(quantile(x$adjusted.results$p0, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$p0, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$p1, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$p1, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$delta, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$delta, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$ratio, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$ratio, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$OR, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$adjusted.results$OR, probs = 1-(1-ci.level)/2, na.rm = TRUE)
      ), ncol=2, byrow=TRUE)
      if (length(x$adjusted.results$p0) == 1) {
        ci_vals <- matrix(rep(NA,10), ncol=2, byrow=TRUE)
      }
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
    tmp <- matrix(c(mean(x$unadjusted.results$p0, na.rm=TRUE), sd(x$unadjusted.results$p0, na.rm=TRUE), mean(x$unadjusted.results$p0, na.rm=TRUE)/sd(x$unadjusted.results$p0, na.rm=TRUE), NA), nrow=1)
    colnames(tmp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp) <- "P0"
    res <- tmp
    
    tmp <- matrix(c(mean(x$unadjusted.results$p1, na.rm=TRUE), sd(x$unadjusted.results$p1, na.rm=TRUE), mean(x$unadjusted.results$p1, na.rm=TRUE)/sd(x$unadjusted.results$p1, na.rm=TRUE), NA), nrow=1)
    rownames(tmp) <- "P1"
    res <- rbind(res,tmp)
    
    tmp <- matrix(c(mean(x$unadjusted.results$delta, na.rm=TRUE), sd(x$unadjusted.results$delta, na.rm=TRUE), mean(x$unadjusted.results$delta, na.rm=TRUE)/sd(x$unadjusted.results$delta, na.rm=TRUE),
                    ifelse(mean(x$unadjusted.results$delta, na.rm=TRUE)/sd(x$unadjusted.results$delta, na.rm=TRUE)<0,
                           2*pnorm(mean(x$unadjusted.results$delta, na.rm=TRUE)/sd(x$unadjusted.results$delta, na.rm=TRUE)),
                           2*(1-pnorm(mean(x$unadjusted.results$delta, na.rm=TRUE)/sd(x$unadjusted.results$delta, na.rm=TRUE))))), nrow=1)
    rownames(tmp) <- "P1-P0"
    res <- rbind(res,tmp)
    
    tmp <- matrix(c(mean(x$unadjusted.results$ratio, na.rm=TRUE), sd(x$unadjusted.results$ratio, na.rm=TRUE), mean(x$unadjusted.results$ratio, na.rm=TRUE)/sd(x$unadjusted.results$ratio, na.rm=TRUE),
                    ifelse(mean(x$unadjusted.results$ratio, na.rm=TRUE)/sd(x$unadjusted.results$ratio, na.rm=TRUE)<0,
                           2*pnorm(mean(x$unadjusted.results$ratio, na.rm=TRUE)/sd(x$unadjusted.results$ratio, na.rm=TRUE)),
                           2*(1-pnorm(mean(x$unadjusted.results$ratio, na.rm=TRUE)/sd(x$unadjusted.results$ratio, na.rm=TRUE))))), nrow=1)
    rownames(tmp) <- "P1/P0"
    res <- rbind(res,tmp)
    
    tmp <- matrix(c(mean(x$unadjusted.results$OR, na.rm=TRUE), sd(x$unadjusted.results$OR, na.rm=TRUE), mean(x$unadjusted.results$OR, na.rm=TRUE)/sd(x$unadjusted.results$OR, na.rm=TRUE),
                    ifelse(mean(x$unadjusted.results$OR, na.rm=TRUE)/sd(x$unadjusted.results$OR, na.rm=TRUE)<0,
                           2*pnorm(mean(x$unadjusted.results$OR, na.rm=TRUE)/sd(x$unadjusted.results$OR, na.rm=TRUE)),
                           2*(1-pnorm(mean(x$unadjusted.results$OR, na.rm=TRUE)/sd(x$unadjusted.results$OR, na.rm=TRUE))))), nrow=1)
    rownames(tmp) <- "OR"
    res <- rbind(res,tmp)
    
    
    if (!is.null(ci.type)) {
      if (ci.type == "norm") {
        ci_vals <- matrix(c(mean(x$unadjusted.results$p0, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$p0, na.rm=TRUE),
                            mean(x$unadjusted.results$p0, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$p0, na.rm=TRUE),
                            mean(x$unadjusted.results$p1, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$p1, na.rm=TRUE),
                            mean(x$unadjusted.results$p1, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$p1, na.rm=TRUE),
                            mean(x$unadjusted.results$delta, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$delta, na.rm=TRUE),
                            mean(x$unadjusted.results$delta, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$delta, na.rm=TRUE),
                            mean(x$unadjusted.results$ratio, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$ratio, na.rm=TRUE),
                            mean(x$unadjusted.results$ratio, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$ratio, na.rm=TRUE),
                            mean(x$unadjusted.results$OR, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$OR, na.rm=TRUE),
                            mean(x$unadjusted.results$OR, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$unadjusted.results$OR, na.rm=TRUE)
        ), ncol=2, byrow=TRUE)
      }
      if (ci.type == "perc") {
        ci_vals <- matrix(c(quantile(x$unadjusted.results$p0, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$p0, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$p1, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$p1, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$delta, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$delta, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$ratio, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$ratio, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$OR, probs = (1-ci.level)/2, na.rm = TRUE),
                            quantile(x$unadjusted.results$OR, probs = 1-(1-ci.level)/2, na.rm = TRUE)
        ), ncol=2, byrow=TRUE)
      }
      tmp <- cbind(res[,1], `Lower CI` = ci_vals[,1], `Upper CI` = ci_vals[,2])
      colnames(tmp)[1] <- "Estimate"
      res <- cbind(tmp, res[,2:4])
    }
    
    printCoefmat(res, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  }
  
  
  cat("\n")
  if (!is.na(object$nevent)) {
    cat(paste0("n= ",x$n,", number of events= ",x$nevent))
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