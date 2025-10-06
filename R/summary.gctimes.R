summary.gctimes <- function (object, digits=4, ci.type=NULL, ci.level=0.95, unadjusted=TRUE, ...)
{
  if (!is.null(ci.type)){
    if (length(ci.type) != 1) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
    if (!(ci.type %in% c("norm","perc"))) {stop("ci.type should be any of the values c(\"norm\",\"perc\")")}
  }
  x <- object
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
  cat("G-computation : \n")
  res_GC <- NULL
  
  mean_AHR <- mean(x$AHR, na.rm=TRUE)
  sd_AHR <- sd(x$AHR, na.rm=TRUE)
  mean_RMST0 <- mean(x$RMST0, na.rm=TRUE)
  sd_RMST0 <- sd(x$RMST0, na.rm=TRUE)
  mean_RMST1 <- mean(x$RMST1, na.rm=TRUE)
  sd_RMST1 <- sd(x$RMST1, na.rm=TRUE)
  mean_deltaRMST <- mean(x$deltaRMST, na.rm=TRUE)
  sd_deltaRMST <- sd(x$deltaRMST, na.rm=TRUE)
  mean_surv0 <- mean(x$surv0, na.rm=TRUE)
  sd_surv0 <- sd(x$surv0, na.rm=TRUE)
  mean_surv1 <- mean(x$surv1, na.rm=TRUE)
  sd_surv1 <- sd(x$surv1, na.rm=TRUE)
  mean_deltasurv <- mean(x$deltasurv, na.rm=TRUE)
  sd_deltasurv <- sd(x$deltasurv, na.rm=TRUE)
  

  if (is.na(sd_AHR) || sd_AHR == 0) {
    z_AHR <- NA
    p_AHR <- NA
  } else {
    z_AHR <- mean_AHR / sd_AHR
    p_AHR <- ifelse(z_AHR < 0, 2 * pnorm(z_AHR), 2 * (1 - pnorm(z_AHR)))
  }
  tmp_AHR <- matrix(c(mean_AHR, sd_AHR, z_AHR, p_AHR), nrow=1)
  colnames(tmp_AHR) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp_AHR) <- "AHR"
  res_GC <- rbind(res_GC, tmp_AHR)
  
  tmp_RMST0 <- matrix(c(mean_RMST0, sd_RMST0, NA, NA), nrow=1)
  colnames(tmp_RMST0) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp_RMST0) <- paste0("RMST0 (to ", x$pro.time,")")
  res_GC <- rbind(res_GC, tmp_RMST0)
  
  tmp_RMST1 <- matrix(c(mean_RMST1, sd_RMST1, NA, NA), nrow=1)
  colnames(tmp_RMST1) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp_RMST1) <- paste0("RMST1 (to ", x$pro.time,")")
  res_GC <- rbind(res_GC, tmp_RMST1)
  
  if (is.na(sd_deltaRMST) || sd_deltaRMST == 0) {
    z_deltaRMST <- NA
    p_deltaRMST <- NA
  } else {
    z_deltaRMST <- mean_deltaRMST / sd_deltaRMST
    p_deltaRMST <- ifelse(z_deltaRMST < 0, 2 * pnorm(z_deltaRMST), 2 * (1 - pnorm(z_deltaRMST)))
  }
  tmp_deltaRMST <- matrix(c(mean_deltaRMST, sd_deltaRMST, z_deltaRMST, p_deltaRMST), nrow=1)
  colnames(tmp_deltaRMST) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp_deltaRMST) <- "RMST1-RMST0"
  res_GC <- rbind(res_GC, tmp_deltaRMST)
  
  tmp_surv0 <- matrix(c(mean_surv0, sd_surv0, NA, NA), nrow=1)
  colnames(tmp_surv0) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp_surv0) <- paste0("S0 (at ", x$pro.time,")")
  res_GC <- rbind(res_GC, tmp_surv0)
  
  tmp_surv1 <- matrix(c(mean_surv1, sd_surv1, NA, NA), nrow=1)
  colnames(tmp_surv1) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp_surv1) <- paste0("S1 (at ", x$pro.time,")")
  res_GC <- rbind(res_GC, tmp_surv1)
  
  if (is.na(sd_deltasurv) || sd_deltasurv == 0) {
    z_deltasurv <- NA
    p_deltasurv <- NA
  } else {
    z_deltasurv <- mean_deltasurv / sd_deltasurv
    p_deltasurv <- ifelse(z_deltasurv < 0, 2 * pnorm(z_deltasurv), 2 * (1 - pnorm(z_deltasurv)))
  }
  tmp_deltasurv <- matrix(c(mean_deltasurv, sd_deltasurv, z_deltasurv, p_deltasurv), nrow=1)
  colnames(tmp_deltasurv) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(tmp_deltasurv) <- paste0("S1 - S0 (at ", x$pro.time,")")
  res_GC <- rbind(res_GC, tmp_deltasurv)
  
  if (!is.null(ci.type)) {
    if (ci.type == "norm") {
      ci_vals_GC <- matrix(c(
        mean_AHR-qnorm(1-(1-ci.level)/2,0,1)*sd_AHR,
        mean_AHR+qnorm(1-(1-ci.level)/2,0,1)*sd_AHR,
        mean_RMST0-qnorm(1-(1-ci.level)/2,0,1)*sd_RMST0,
        mean_RMST0+qnorm(1-(1-ci.level)/2,0,1)*sd_RMST0,
        mean_RMST1-qnorm(1-(1-ci.level)/2,0,1)*sd_RMST1,
        mean_RMST1+qnorm(1-(1-ci.level)/2,0,1)*sd_RMST1,
        mean_deltaRMST-qnorm(1-(1-ci.level)/2,0,1)*sd_deltaRMST,
        mean_deltaRMST+qnorm(1-(1-ci.level)/2,0,1)*sd_deltaRMST,
        mean_surv0-qnorm(1-(1-ci.level)/2,0,1)*sd_surv0,
        mean_surv0+qnorm(1-(1-ci.level)/2,0,1)*sd_surv0,
        mean_surv1-qnorm(1-(1-ci.level)/2,0,1)*sd_surv1,
        mean_surv1+qnorm(1-(1-ci.level)/2,0,1)*sd_surv1,
        mean_deltasurv-qnorm(1-(1-ci.level)/2,0,1)*sd_deltasurv,
        mean_deltasurv+qnorm(1-(1-ci.level)/2,0,1)*sd_deltasurv
      ), ncol=2, byrow=TRUE)
    }
    if (ci.type == "perc") {
      ci_vals_GC <- matrix(c(
        quantile(x$AHR, probs = (1-ci.level)/2, na.rm = TRUE),
        quantile(x$AHR, probs = 1-(1-ci.level)/2, na.rm = TRUE),
        quantile(x$RMST0, probs = (1-ci.level)/2, na.rm = TRUE),
        quantile(x$RMST0, probs = 1-(1-ci.level)/2, na.rm = TRUE),
        quantile(x$RMST1, probs = (1-ci.level)/2, na.rm = TRUE),
        quantile(x$RMST1, probs = 1-(1-ci.level)/2, na.rm = TRUE),
        quantile(x$deltaRMST, probs = (1-ci.level)/2, na.rm = TRUE),
        quantile(x$deltaRMST, probs = 1-(1-ci.level)/2, na.rm = TRUE),
        quantile(x$surv0, probs = (1-ci.level)/2, na.rm = TRUE),
        quantile(x$surv0, probs = 1-(1-ci.level)/2, na.rm = TRUE),
        quantile(x$surv1, probs = (1-ci.level)/2, na.rm = TRUE),
        quantile(x$surv1, probs = 1-(1-ci.level)/2, na.rm = TRUE),
        quantile(x$deltasurv, probs = (1-ci.level)/2, na.rm = TRUE),
        quantile(x$deltasurv, probs = 1-(1-ci.level)/2, na.rm = TRUE)
      ), ncol=2, byrow=TRUE)
    }
    tmp_GC <- cbind(res_GC[,1], `Lower CI` = ci_vals_GC[,1], `Upper CI` = ci_vals_GC[,2])
    colnames(tmp_GC)[1] <- "Estimate"
    res_GC <- cbind(tmp_GC, res_GC[,2:4])
  }
  printCoefmat(res_GC, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  cat("\n")
  res_GC_return <- res_GC
  
  
  if (!is.null(object$newdata)) {unadjusted = FALSE}
  if (unadjusted == TRUE) {
    cat("Unadjusted : \n")
    res_unadj <- NULL
    
    mean_AHR_unadj <- mean(x$AHR.unadj, na.rm=TRUE)
    sd_AHR_unadj <- sd(x$AHR.unadj, na.rm=TRUE)
    mean_RMST0_unadj <- mean(x$RMST0.unadj, na.rm=TRUE)
    sd_RMST0_unadj <- sd(x$RMST0.unadj, na.rm=TRUE)
    mean_RMST1_unadj <- mean(x$RMST1.unadj, na.rm=TRUE)
    sd_RMST1_unadj <- sd(x$RMST1.unadj, na.rm=TRUE)
    mean_deltaRMST_unadj <- mean(x$deltaRMST.unadj, na.rm=TRUE)
    sd_deltaRMST_unadj <- sd(x$deltaRMST.unadj, na.rm=TRUE)
    mean_surv0_unadj <- mean(x$surv0.unadj, na.rm=TRUE)
    sd_surv0_unadj <- sd(x$surv0.unadj, na.rm=TRUE)
    mean_surv1_unadj <- mean(x$surv1.unadj, na.rm=TRUE)
    sd_surv1_unadj <- sd(x$surv1.unadj, na.rm=TRUE)
    mean_deltasurv_unadj <- mean(x$deltasurv.unadj, na.rm=TRUE)
    sd_deltasurv_unadj <- sd(x$deltasurv.unadj, na.rm=TRUE)
    
    
    if (is.na(sd_AHR_unadj) || sd_AHR_unadj == 0) {
      z_AHR_unadj <- NA
      p_AHR_unadj <- NA
    } else {
      z_AHR_unadj <- mean_AHR_unadj / sd_AHR_unadj
      p_AHR_unadj <- ifelse(z_AHR_unadj < 0, 2 * pnorm(z_AHR_unadj), 2 * (1 - pnorm(z_AHR_unadj)))
    }
    tmp_AHR_unadj <- matrix(c(mean_AHR_unadj, sd_AHR_unadj, z_AHR_unadj, p_AHR_unadj), nrow=1)
    colnames(tmp_AHR_unadj) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp_AHR_unadj) <- "AHR"
    res_unadj <- rbind(res_unadj, tmp_AHR_unadj)
    
    tmp_RMST0_unadj <- matrix(c(mean_RMST0_unadj, sd_RMST0_unadj, NA, NA), nrow=1)
    colnames(tmp_RMST0_unadj) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp_RMST0_unadj) <- paste0("RMST0 (to ", x$pro.time,")")
    res_unadj <- rbind(res_unadj, tmp_RMST0_unadj)
    
    tmp_RMST1_unadj <- matrix(c(mean_RMST1_unadj, sd_RMST1_unadj, NA, NA), nrow=1)
    colnames(tmp_RMST1_unadj) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp_RMST1_unadj) <- paste0("RMST1 (to ", x$pro.time,")")
    res_unadj <- rbind(res_unadj, tmp_RMST1_unadj)
    
    if (is.na(sd_deltaRMST_unadj) || sd_deltaRMST_unadj == 0) {
      z_deltaRMST_unadj <- NA
      p_deltaRMST_unadj <- NA
    } else {
      z_deltaRMST_unadj <- mean_deltaRMST_unadj / sd_deltaRMST_unadj
      p_deltaRMST_unadj <- ifelse(z_deltaRMST_unadj < 0, 2 * pnorm(z_deltaRMST_unadj), 2 * (1 - pnorm(z_deltaRMST_unadj)))
    }
    tmp_deltaRMST_unadj <- matrix(c(mean_deltaRMST_unadj, sd_deltaRMST_unadj, z_deltaRMST_unadj, p_deltaRMST_unadj), nrow=1)
    colnames(tmp_deltaRMST_unadj) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp_deltaRMST_unadj) <- "RMST1-RMST0"
    res_unadj <- rbind(res_unadj, tmp_deltaRMST_unadj)
    
    tmp_surv0_unadj <- matrix(c(mean_surv0_unadj, sd_surv0_unadj, NA, NA), nrow=1)
    colnames(tmp_surv0_unadj) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp_surv0_unadj) <- paste0("S0 (at ", x$pro.time,")")
    res_unadj <- rbind(res_unadj, tmp_surv0_unadj)
    
    tmp_surv1_unadj <- matrix(c(mean_surv1_unadj, sd_surv1_unadj, NA, NA), nrow=1)
    colnames(tmp_surv1_unadj) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp_surv1_unadj) <- paste0("S1 (at ", x$pro.time,")")
    res_unadj <- rbind(res_unadj, tmp_surv1_unadj)
    
    if (is.na(sd_deltasurv_unadj) || sd_deltasurv_unadj == 0) {
      z_deltasurv_unadj <- NA
      p_deltasurv_unadj <- NA
    } else {
      z_deltasurv_unadj <- mean_deltasurv_unadj / sd_deltasurv_unadj
      p_deltasurv_unadj <- ifelse(z_deltasurv_unadj < 0, 2 * pnorm(z_deltasurv_unadj), 2 * (1 - pnorm(z_deltasurv_unadj)))
    }
    tmp_deltasurv_unadj <- matrix(c(mean_deltasurv_unadj, sd_deltasurv_unadj, z_deltasurv_unadj, p_deltasurv_unadj), nrow=1)
    colnames(tmp_deltasurv_unadj) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
    rownames(tmp_deltasurv_unadj) <- paste0("S1 - S0 (at ", x$pro.time,")")
    res_unadj <- rbind(res_unadj, tmp_deltasurv_unadj)
    
    if (!is.null(ci.type)) {
      if (ci.type == "norm") {
        ci_vals_unadj <- matrix(c(
          mean_AHR_unadj-qnorm(1-(1-ci.level)/2,0,1)*sd_AHR_unadj,
          mean_AHR_unadj+qnorm(1-(1-ci.level)/2,0,1)*sd_AHR_unadj,
          mean_RMST0_unadj-qnorm(1-(1-ci.level)/2,0,1)*sd_RMST0_unadj,
          mean_RMST0_unadj+qnorm(1-(1-ci.level)/2,0,1)*sd_RMST0_unadj,
          mean_RMST1_unadj-qnorm(1-(1-ci.level)/2,0,1)*sd_RMST1_unadj,
          mean_RMST1_unadj+qnorm(1-(1-ci.level)/2,0,1)*sd_RMST1_unadj,
          mean_deltaRMST_unadj-qnorm(1-(1-ci.level)/2,0,1)*sd_deltaRMST_unadj,
          mean_deltaRMST_unadj+qnorm(1-(1-ci.level)/2,0,1)*sd_deltaRMST_unadj,
          mean_surv0_unadj-qnorm(1-(1-ci.level)/2,0,1)*sd_surv0_unadj,
          mean_surv0_unadj+qnorm(1-(1-ci.level)/2,0,1)*sd_surv0_unadj,
          mean_surv1_unadj-qnorm(1-(1-ci.level)/2,0,1)*sd_surv1_unadj,
          mean_surv1_unadj+qnorm(1-(1-ci.level)/2,0,1)*sd_surv1_unadj,
          mean_deltasurv_unadj-qnorm(1-(1-ci.level)/2,0,1)*sd_deltasurv_unadj,
          mean_deltasurv_unadj+qnorm(1-(1-ci.level)/2,0,1)*sd_deltasurv_unadj
        ), ncol=2, byrow=TRUE)
      }
      if (ci.type == "perc") {
        ci_vals_unadj <- matrix(c(
          quantile(x$AHR.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
          quantile(x$AHR.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
          quantile(x$RMST0.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
          quantile(x$RMST0.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
          quantile(x$RMST1.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
          quantile(x$RMST1.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
          quantile(x$deltaRMST.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
          quantile(x$deltaRMST.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
          quantile(x$surv0.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
          quantile(x$surv0.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
          quantile(x$surv1.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
          quantile(x$surv1.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE),
          quantile(x$deltasurv.unadj, probs = (1-ci.level)/2, na.rm = TRUE),
          quantile(x$deltasurv.unadj, probs = 1-(1-ci.level)/2, na.rm = TRUE)
        ), ncol=2, byrow=TRUE)
      }
      tmp_unadj <- cbind(res_unadj[,1], `Lower CI` = ci_vals_unadj[,1], `Upper CI` = ci_vals_unadj[,2])
      colnames(tmp_unadj)[1] <- "Estimate"
      res_unadj <- cbind(tmp_unadj, res_unadj[,2:4])
    }
    printCoefmat(res_unadj, digits = digits, P.values = TRUE, has.Pvalue = TRUE, na.print = "", ...)
  }
  
  cat("\n")
  if (!is.na(object$nevent)) {
    cat(paste0("n= ",x$n,", number of events= ",x$nevent))
  } else {
    cat(paste0("n= ",x$n))
  }
  cat("\n")
  if(x$missing==1) { cat(x$missing, " observation deleted due to missingness", sep=""); cat("\n") }
  if(x$missing >1) { cat(x$missing, " observations deleted due to missingness", sep=""); cat("\n") }
  
  if (unadjusted == TRUE) {
    invisible(list(
      GC = as.data.frame(res_GC_return),
      unadjusted = as.data.frame(res_unadj)))
  } else {
    invisible(list(
      GC = as.data.frame(res_GC_return)))
  }
}