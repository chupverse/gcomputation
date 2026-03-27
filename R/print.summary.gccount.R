

print.summary.gccount <- function (x, ...)
{
  digits <- x$digits
  if (x$model %in% c("lasso","ridge","elasticnet")) {cat("model: ",x$model,", tuning parameters: ", sep = "")
    if (x$model == "elasticnet") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep=" ")}
    if (x$model == "lasso") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep=" ")}
    if (x$model == "ridge") {
      cat("lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep=" ")}
    cat("\nfamily: poisson(link = 'log')")
    cat("\nCall:", "\n", sep = "")
    dput(x$formula)}
  if (x$model %in% c("all","aic","bic")) {cat(x$model," model \nCall:", "\n", sep = "")
    if (is.null(x$m)) {dput(x$tuning.parameters)} else {dput(x$tuning.parameters[[1]])}
  }
  cat("\n")
  
  cat("G-computation: \n")
  printCoefmat(x$adjusted, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  cat("\n")
  
  if (x$unadjusted.flag == TRUE) {
    cat("Unadjusted: \n")
    printCoefmat(x$unadjusted, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  }
  
  cat("\n")
  cat(paste0("n= ",x$n))
  cat("\n")
  if (!is.null(x$nimput)) {
    if (x$nimput == 1) { cat(x$nimput, " observation imputed", sep=""); cat("\n") }
    if (x$nimput > 1) { cat(x$nimput, " observations imputed", sep=""); cat("\n") }
  }
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep=""); cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep=""); cat("\n") }
  
  invisible(x)
}