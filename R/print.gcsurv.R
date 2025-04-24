print.gcsurv <- function (x, digits=4, ...)
{
  cat("Call:", "\n", sep = "")
  dput(x$formula)
  cat("\n")
  
  tmp <- matrix(c(x$RMST$mean.deltaRMST, x$RMST$se.deltaRMST,
                  x$RMST$ci.low.asympt.deltaRMST, x$RMST$ci.upp.asympt.deltaRMST,
                  x$RMST$mean.deltaRMST/x$RMST$se.deltaRMST, x$RMST$p.value.deltaRMST),nrow=1)
  colnames(tmp) <- c("deltaRMST","se(deltaRMST)","lower .95(asympt.)", "upper .95(asympt.)","z", "p")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, ...)
  cat("\n")
  
  cat(paste0("n= ",nrow(x$data),", number of events= ",sum(x$data[,x$outcomes[[2]]])))
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep="") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep="") }
}