print.gcbinary <- function (x, digits=4, ...)
{
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
  
  cat("Estimates : \n")
  res <- matrix(c(mean(x$p0, na.rm=TRUE),
      mean(x$p1, na.rm=TRUE),
      mean(x$delta, na.rm=TRUE),
      mean(x$ratio, na.rm=TRUE),
      mean(x$OR, na.rm=TRUE)),
    nrow = 1)
  colnames(res) <- c("P0", "P1", "P1-P0", "P1/P0", "OR")
  rownames(res) <-  ""
  
  printCoefmat(res, digits = digits, ..., P.values = FALSE, has.Pvalue = FALSE, na.print = "")
  
  
  cat("\n")
  
  if (!is.na(x$nevent)) {
    cat(paste0("n= ",x$n,", number of events= ",x$nevent))
  } else {
    cat(paste0("n= ",x$n))
  }
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep=""); cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep=""); cat("\n") }
}