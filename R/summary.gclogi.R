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
  
  tmp <- matrix(c(x$coefficients$mean.delta, x$coefficients$se.delta,
                  x$coefficients$ci.low.asympt.delta, x$coefficients$ci.upp.asympt.delta,
                  x$coefficients$ci.low.nonpara.delta, x$coefficients$ci.upp.nonpara.delta, x$coefficients$mean.delta/x$coefficients$se.delta, x$coefficients$p.value.delta),nrow=1)
  colnames(tmp) <- c("delta","se(delta)","lower .95(asympt.)", "upper .95(asympt.)","lower .95(non para.)", "upper .95(non para.)","z", "Pr(>|z|)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$coefficients$mean.p0, x$coefficients$se.p0,
                  x$coefficients$ci.low.asympt.p0, x$coefficients$ci.upp.asympt.p0,
                  x$coefficients$ci.low.nonpara.p0, x$coefficients$ci.upp.nonpara.p0),nrow=1)
  colnames(tmp) <- c("p0","se(p0)","lower .95(asympt.)", "upper .95(asympt.)","lower .95(non para.)", "upper .95(non para.)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, ...)
  cat("\n")
  
  tmp <- matrix(c(x$coefficients$mean.p1, x$coefficients$se.p1,
                  x$coefficients$ci.low.asympt.p1, x$coefficients$ci.upp.asympt.p1,
                  x$coefficients$ci.low.nonpara.p1, x$coefficients$ci.upp.nonpara.p1),nrow=1)
  colnames(tmp) <- c("p1","se(p1)","lower .95(asympt.)", "upper .95(asympt.)","lower .95(non para.)", "upper .95(non para.)")
  rownames(tmp) = ""
  printCoefmat(tmp, digits = digits, P.values = TRUE, 
               has.Pvalue = TRUE, ...)
  cat("\n")
  
  if (x$method == "elasticnet") {
    cat("Tuning parameters: lambda= ",round(x$tuning.parameters$lambda,digits=digits), " alpha= ",round(x$tuning.parameters$alpha,digits=digits),sep="")}
  if (x$method == "ridge" | x$method == "ridge" ) {
    cat("Tuning parameters: lambda= ",round(x$tuning.parameters$lambda,digits=digits),sep="")}
  
  
}