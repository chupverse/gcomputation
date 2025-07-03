summary.gcsurv <- function (object, digits=4, ...)
{
  x <- object
  cat("Call:", "\n", sep = "")
  dput(x$formula)
  cat("\n")
  cat(paste0("n= ",nrow(x$data),", number of events= ",sum(x$data[,x$outcome[[2]]])))
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness \n", sep="") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness \n", sep="") }
  cat("\n")
  
  tmp <- matrix(c(x$RMST$mean.deltaRMST, x$RMST$se.deltaRMST,
                  x$RMST$ci.low.asympt.deltaRMST, x$RMST$ci.upp.asympt.deltaRMST,
                  x$RMST$ci.low.nonpara.deltaRMST, x$RMST$ci.upp.nonpara.deltaRMST, x$RMST$mean.deltaRMST/x$RMST$se.deltaRMST, x$RMST$p.value.deltaRMST),nrow=1)
  colnames(tmp) <- c("deltaRMST","se(deltaRMST)","lower .95(asympt.)", "upper .95(asympt.)","lower .95(non para.)", "upper .95(non para.)","z", "Pr(>|z|)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$RMST$mean.RMST0, x$RMST$se.RMST0,
                  x$RMST$ci.low.asympt.RMST0, x$RMST$ci.upp.asympt.RMST0,
                  x$RMST$ci.low.nonpara.RMST0, x$RMST$ci.upp.nonpara.RMST0),nrow=1)
  colnames(tmp) <- c("RMST0","se(RMST0)","lower .95(asympt.)", "upper .95(asympt.)","lower .95(non para.)", "upper .95(non para.)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$RMST$mean.RMST1, x$RMST$se.RMST1,
                  x$RMST$ci.low.asympt.RMST1, x$RMST$ci.upp.asympt.RMST1,
                  x$RMST$ci.low.nonpara.RMST1, x$RMST$ci.upp.nonpara.RMST1),nrow=1)
  colnames(tmp) <- c("RMST1","se(RMST1)","lower .95(asympt.)", "upper .95(asympt.)","lower .95(non para.)", "upper .95(non para.)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, ...)
  cat("\n")
  
  if (x$method == "elasticnet") {
    cat("Tuning parameters: lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep="")}
  if (x$method == "ridge" | x$method == "ridge" ) {
    cat("Tuning parameters: lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep="")}
  
  
}