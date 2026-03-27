print.gcbinary <- function (x, digits=4, ...)
{
  if (x$model %in% c("lasso","ridge","elasticnet")) {
    cat("model: ",x$model, sep = "")
    if (is.null(x$m)) {
      cat(", tuning parameters: ", sep = "")
      if (x$model == "elasticnet") {
        cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep=" ")}
      if (x$model == "lasso") {
        cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep=" ")}
      if (x$model == "ridge") {
        cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep=" ")}
    }
    cat("\nfamily: binomial(link = 'logit')")
    cat("\nCall:", "\n", sep = "")
    dput(x$formula)
  }
  if (x$model %in% c("all","aic","bic")) {
    cat(x$model," model \nfamily: binomial(link = 'logit')\nCall:", "\n", sep = "")
    if (is.null(x$m)) {dput(x$tuning.parameters)} else {dput(x$tuning.parameters[[1]])}
  }
  cat("\n")
  
  cat("Estimates: \n")
  res <- matrix(c(mean(x$adjusted.results$p0, na.rm=TRUE),
                  mean(x$adjusted.results$p1, na.rm=TRUE),
                  mean(x$adjusted.results$delta, na.rm=TRUE),
                  mean(x$adjusted.results$ratio, na.rm=TRUE),
                  mean(x$adjusted.results$OR, na.rm=TRUE)),
                nrow = 1)
  colnames(res) <- c("p0", "p1", "p1-p0", "p1/p0", "OR")
  rownames(res) <-  ""
  
  printCoefmat(res, digits = digits, ..., P.values = FALSE, has.Pvalue = FALSE, na.print = "")
  
  
  cat("\n")
  
  if (!is.na(x$nevent)) {
    cat(paste0("n= ",x$n,", number of events= ",x$nevent))
  } else {
    cat(paste0("n= ",x$n))
  }
  cat("\n")
  if (!is.null(x$nimput)) {
    if (x$nimput == 1) { cat(x$nimput, " observation imputed", sep=""); cat("\n") }
    if (x$nimput > 1) { cat(x$nimput, " observations imputed", sep=""); cat("\n") }
  }
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep=""); cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep=""); cat("\n") }
}