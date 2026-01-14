print.gctimes <- function (x, digits=4, ...)
{
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
  
  cat("Estimates : \n")
  res <- matrix(c( mean(x$adjusted.results$AHR, na.rm=TRUE),
                   mean(x$adjusted.results$RMST0, na.rm=TRUE),
                   mean(x$adjusted.results$RMST1, na.rm=TRUE),
                   mean(x$adjusted.results$deltaRMST, na.rm=TRUE),
                   mean(x$adjusted.results$s0, na.rm=TRUE),
                   mean(x$adjusted.results$s1, na.rm=TRUE),
                   mean(x$adjusted.results$delta, na.rm=TRUE)
  ), nrow = 1)
  
  colnames(res) <- c("AHR", 
                     paste0("RMST0(to ", x$pro.time,")"),
                     paste0("RMST1(to ", x$pro.time,")"),
                     "RMST1-RMST0",
                     paste0("S0(at ", x$pro.time,")"),
                     paste0("S1(at ", x$pro.time,")"),
                     paste0("S1 - S0(at ", x$pro.time,")"))
  rownames(res) <- ""
  
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
