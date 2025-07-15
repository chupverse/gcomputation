# Pq il y a "if (x$method == "elasticnet")" mais pas la même condition avec "lasso"

summary.gclogi <- function (object, digits=4, ...)
{
  x <- object
  cat("Call:", "\n", sep = "")
  dput(x$formula)
  cat("\n")
  cat(paste0("n= ",nrow(x$data),", number of events= ",sum(x$data[,x$outcome])))
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness \n", sep="") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness \n", sep="") }
  cat("\n")
  
  tmp <- matrix(c(x$delta$estim, x$delta$se, x$delta$ci.asympt[1], x$delta$ci.asympt[2],
                  x$delta$ci.nonpara[1], x$delta$ci.nonpara[2], x$delta$estim/x$delta$se, x$delta$p.value), nrow=1)
  colnames(tmp) <- c("delta","se(delta)","lower .95(asympt.)", "upper .95(asympt.)", "lower .95(non para.)", "upper .95(non para.)", "z", "Pr(>|z|)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, has.Pvalue = TRUE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$p0$estim, x$p0$se, x$p0$ci.asympt[1], x$p0$ci.asympt[2], x$p0$ci.nonpara[1], x$p0$ci.nonpara[2]), nrow=1)
  colnames(tmp) <- c("p0", "se(p0)", "lower .95(asympt.)", "upper .95(asympt.)", "lower .95(non para.)", "upper .95(non para.)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, has.Pvalue = TRUE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$p1$estim, x$p1$se, x$p1$ci.asympt[1], x$p1$ci.asympt[2], x$p1$ci.nonpara[1], x$p1$ci.nonpara[2]), nrow=1)
  colnames(tmp) <- c("p1", "se(p1)", "lower .95(asympt.)", "upper .95(asympt.)", "lower .95(non para.)", "upper .95(non para.)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, has.Pvalue = TRUE, ...)
  cat("\n")
  
  if (x$method == "elasticnet") {
    cat("Tuning parameters: lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep="")}
  if (x$method == "ridge" | x$method == "ridge" ) {
    cat("Tuning parameters: lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep="")}
}