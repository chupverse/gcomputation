summary.gclogi <- function (object, digits=4, ci.type=NULL, ci.level=0.95, ...)
{
  if (!is.null(ci.type)){
    if (length(ci.type) != 1) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
    if (!(ci.type %in% c("norm","perc"))) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
  }

  x <- object
  if (x$method %in% c("lasso","ridge","elasticnet")) {cat("Method: ",x$method,", tuning parameters: ", sep = "")
    if (x$method == "elasticnet") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep="")}
    if (x$method == "lasso") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep="")}
    if (x$method == "ridge") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep="")}
    cat("\nCall:", "\n", sep = "")
    dput(x$formula)}
  if (x$method %in% c("all","aic","bic")) {cat(x$method," method \nCall:", "\n", sep = "")
    dput(x$tuning.parameters)}
  cat("\n")

  
  cat("G-computation : \n")
  tmp <- matrix(c(mean(x$p0, na.rm=TRUE), sd(x$p0, na.rm=TRUE), mean(x$p0, na.rm=TRUE)/sd(x$p0, na.rm=TRUE), NA), nrow=1)
  colnames(tmp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp) <- "P0"
  res <- tmp
  
  tmp <- matrix(c(mean(x$p1, na.rm=TRUE), sd(x$p1, na.rm=TRUE), mean(x$p1, na.rm=TRUE)/sd(x$p1, na.rm=TRUE), NA), nrow=1)
  rownames(tmp) <- "P1"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$delta, na.rm=TRUE), sd(x$delta, na.rm=TRUE), mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE),
                  ifelse(mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE)<0,
                         2*pnorm(mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$delta, na.rm=TRUE)/sd(x$delta, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "Difference (1-0)"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$ratio, na.rm=TRUE), sd(x$ratio, na.rm=TRUE), mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE),
                  ifelse(mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE)<0,
                         2*pnorm(mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$ratio, na.rm=TRUE)/sd(x$ratio, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "Ratio (1/0)"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$OR, na.rm=TRUE), sd(x$OR, na.rm=TRUE), mean(x$OR, na.rm=TRUE)/sd(x$OR, na.rm=TRUE),
                  ifelse(mean(x$OR, na.rm=TRUE)/sd(x$OR, na.rm=TRUE)<0,
                         2*pnorm(mean(x$OR, na.rm=TRUE)/sd(x$OR, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$OR, na.rm=TRUE)/sd(x$OR, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "OR"
  res <- rbind(res,tmp)
  
  
  if (!is.null(ci.type)) {
    if (ci.type == "norm") {
      ci_vals <- matrix(c(mean(x$p0, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p0, na.rm=TRUE),
                          mean(x$p0, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p0, na.rm=TRUE),
                          mean(x$p1, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p1, na.rm=TRUE),
                          mean(x$p1, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p1, na.rm=TRUE),
                          mean(x$delta, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta, na.rm=TRUE),
                          mean(x$delta, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta, na.rm=TRUE),
                          mean(x$ratio, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio, na.rm=TRUE),
                          mean(x$ratio, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio, na.rm=TRUE),
                          mean(x$OR, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$OR, na.rm=TRUE),
                          mean(x$OR, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$OR, na.rm=TRUE)
                          ), ncol=2, byrow=TRUE)
    }
    if (ci.type == "perc") {
      ci_vals <- matrix(c(quantile(x$p0, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$p0, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$p1, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$p1, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$delta, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$delta, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$ratio, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$ratio, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$OR, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$OR, probs = 1-(1-ci.level)/2, na.rm = TRUE)
      ), ncol=2, byrow=TRUE)
    }
    tmp <- cbind(res[,1], `Lower CI` = ci_vals[,1], `Upper CI` = ci_vals[,2])
    colnames(tmp)[1] <- "Estimate"
    res <- cbind(tmp, res[,2:4])
  }
  
  printCoefmat(res, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  cat("\n")
  
  res_GC <- res
  
  cat("Unadjusted : \n")
  tmp <- matrix(c(mean(x$p0.unadj, na.rm=TRUE), sd(x$p0.unadj, na.rm=TRUE), mean(x$p0.unadj, na.rm=TRUE)/sd(x$p0.unadj, na.rm=TRUE), NA), nrow=1)
  colnames(tmp) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp) <- "P0"
  res <- tmp
  
  tmp <- matrix(c(mean(x$p1.unadj, na.rm=TRUE), sd(x$p1.unadj, na.rm=TRUE), mean(x$p1.unadj, na.rm=TRUE)/sd(x$p1.unadj, na.rm=TRUE), NA), nrow=1)
  rownames(tmp) <- "P1"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$delta.unadj, na.rm=TRUE), sd(x$delta.unadj, na.rm=TRUE), mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE),
                  ifelse(mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE)<0,
                         2*pnorm(mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$delta.unadj, na.rm=TRUE)/sd(x$delta.unadj, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "Difference (1-0)"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$ratio.unadj, na.rm=TRUE), sd(x$ratio.unadj, na.rm=TRUE), mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE),
                  ifelse(mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE)<0,
                         2*pnorm(mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$ratio.unadj, na.rm=TRUE)/sd(x$ratio.unadj, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "Ratio (1/0)"
  res <- rbind(res,tmp)
  
  tmp <- matrix(c(mean(x$OR.unadj, na.rm=TRUE), sd(x$OR.unadj, na.rm=TRUE), mean(x$OR.unadj, na.rm=TRUE)/sd(x$OR.unadj, na.rm=TRUE),
                  ifelse(mean(x$OR.unadj, na.rm=TRUE)/sd(x$OR.unadj, na.rm=TRUE)<0,
                         2*pnorm(mean(x$OR.unadj, na.rm=TRUE)/sd(x$OR.unadj, na.rm=TRUE)),
                         2*(1-pnorm(mean(x$OR.unadj, na.rm=TRUE)/sd(x$OR.unadj, na.rm=TRUE))))), nrow=1)
  rownames(tmp) <- "OR"
  res <- rbind(res,tmp)
  
  
  if (!is.null(ci.type)) {
    if (ci.type == "norm") {
      ci_vals <- matrix(c(mean(x$p0.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p0.unadj, na.rm=TRUE),
                          mean(x$p0.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p0.unadj, na.rm=TRUE),
                          mean(x$p1.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p1.unadj, na.rm=TRUE),
                          mean(x$p1.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$p1.unadj, na.rm=TRUE),
                          mean(x$delta.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta.unadj, na.rm=TRUE),
                          mean(x$delta.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$delta.unadj, na.rm=TRUE),
                          mean(x$ratio.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio.unadj, na.rm=TRUE),
                          mean(x$ratio.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$ratio.unadj, na.rm=TRUE),
                          mean(x$OR.unadj, na.rm=TRUE) - qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$OR.unadj, na.rm=TRUE),
                          mean(x$OR.unadj, na.rm=TRUE) + qnorm(1-(1-ci.level)/2, 0, 1)*sd(x$OR.unadj, na.rm=TRUE)
      ), ncol=2, byrow=TRUE)
    }
    if (ci.type == "perc") {
      ci_vals <- matrix(c(quantile(x$p0.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$p0.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$p1.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$p1.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$delta.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$delta.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$ratio.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$ratio.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
                          quantile(x$OR.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
                          quantile(x$OR.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE)
      ), ncol=2, byrow=TRUE)
    }
    tmp <- cbind(res[,1], `Lower CI` = ci_vals[,1], `Upper CI` = ci_vals[,2])
    colnames(tmp)[1] <- "Estimate"
    res <- cbind(tmp, res[,2:4])
  }
  
  printCoefmat(res, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  cat("\n")
  
  
  cat(paste0("n= ",x$n,", number of events= ",x$nevent))
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep=""); cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep=""); cat("\n") }
  
  invisible(list(
    GC = res_GC,
    unadjusted = res
  ))
}