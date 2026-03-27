

print.summary.gctimes <- function (x, ...)
{
  digits <- x$digits
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
    cat("\nProportional hazard regression (Breslow estimator)")
    cat("\nCall:", "\n", sep = "")
    dput(x$formula)
  }
  if (x$model %in% c("all","aic","bic")) {
    cat(x$model," model \nProportional hazard regression (Breslow estimator)\nCall:", "\n", sep = "")
    if (is.null(x$m)) {dput(x$tuning.parameters)} else {dput(x$tuning.parameters[[1]])}
  }
  cat("\n")
  cat("G-computation: \n")
  printCoefmat(as.matrix(x$adjusted), digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  cat("\n")
  
  if (x$unadjusted.flag == TRUE) {
    cat("Unadjusted: \n")
    printCoefmat(as.matrix(x$unadjusted), digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  }
  
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
  
  invisible(x)
}