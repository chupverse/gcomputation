print.gclogi <- function (x, digits=4, ...)
{
  cat("Call:", "\n", sep = "")
  dput(x$formula)
  cat("\n")
  
  tmp <- matrix(c(x$coefficient$mean.delta, x$coefficient$se.delta,
                  x$coefficient$mean.delta/x$coefficient$se.delta, x$coefficient$p.value.delta),nrow=1)
  colnames(tmp) <- c("delta","se(delta)","z", "p")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, signif.stars = FALSE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$mOR$mean.mOR, x$mOR$se.mOR,
                  x$mOR$mean.mOR/x$mOR$se.mOR, x$mOR$p.value.mOR),nrow=1)
  colnames(tmp) <- c("mOR","se(mOR)","z", "p")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, signif.stars = FALSE, ...)
  cat("\n")
  
  cat(paste0("n= ",x$n,", number of events= ",x$nevent))
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep="");cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep="");cat("\n") }
}