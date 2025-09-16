print.gclogi <- function (x, digits=4, ...)
{
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
  
  cat("Estimates : \n")
  res <- matrix(c(mean(x$p0, na.rm=TRUE),
      mean(x$p1, na.rm=TRUE),
      mean(x$delta, na.rm=TRUE),
      mean(x$ratio, na.rm=TRUE),
      mean(x$OR, na.rm=TRUE)),
    nrow = 1)
  colnames(res) <- c("P0", "P1", "Difference(1-0)", "Ratio(1/0)", "OR")
  rownames(res) <-  ""
  
  printCoefmat(res, digits = digits, ..., P.values = FALSE, has.Pvalue = FALSE, na.print = "")
  
  
  cat("\n")
  
  cat(paste0("n= ",x$n,", number of events= ",x$nevent))
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep=""); cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep=""); cat("\n") }
}