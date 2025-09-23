plot.gctimes <- function (x, method="calibration", n.groups=5, pro.time=NULL, smooth=FALSE, ...) {
  if (!(method %in% c("calibration","survival"))) {stop("Method needs to be calibration or survival")}
  if (method == "survival") {
    data = x$data
    times = x$outcome$times
    failures = x$outcome$failures
    results.surv.calibration = x$calibration
    
    km <- survfit(formula(paste0("Surv(",times,",", failures,")~1")), data=data) 
    ## Duplicate rows if more than one event to be able to compare with GC (Basically remove any where no events!)
    event_rows <- km$n.event > 0  
    reptimes <- km$time[event_rows]
    repcumhaz <- km$cumhaz[event_rows]
    repsurv <- km$surv[event_rows]
    
    results.km <- data.frame(time=c(0, reptimes), cumhaz=c(0, repcumhaz), surv=c(1,repsurv) ) #return
    
    

      plot(c(0, results.km$time), c(1, results.km$surv),ylab="Survival Probability", xlab="Time", type="s", col=adjustcolor("blue", alpha.f = 0.6), lwd=2, ...)
    lines(results.surv.calibration$time, results.surv.calibration$surv,
          col=adjustcolor("red", alpha.f = 0.6), lwd=2, lty=1, type="s") #add ifelse type = "l" if new non parametric method added to function
    grid()
    legend("topright", legend=c("Kaplan-Meier Estimate", "Mean of survival predictions"), col=c("blue", "red"), lty=c(1, 1), lwd=2, bty="n")
  }
  
  
  if (method == "calibration") {
    
    
    if (!is.null(x$m)) {
      if (smooth) {
        all_pred <- c()
        all_obs <- c()
        for (i in 1:x$m) {
          if(is.null(pro.time)) pro.time <- median(x$data[[i]][,x$outcome[[1]]])
          time <- x$data[[i]][,x$outcome[[1]]]
          event <- x$data[[i]][,x$outcome[[2]]]
          .lp <- x$calibration[[i]]$lp
          hazard <- x$calibration[[i]]$H0.multi
          T.multi <- x$calibration[[i]]$time
          
          .pred <- exp(matrix(exp(.lp)) %*% t(as.matrix(-1*hazard)))
          pos <- max(which(T.multi <= pro.time))
          .pred <- .pred[,pos]
          
          .grps <- as.numeric(cut(.pred,
                                  breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                                  labels = 1:n.groups))
          
          .est <- sapply(1:n.groups, function(x) mean(.pred[.grps==x]))
          
          .data <- data.frame(time=time, event=event, grps=.grps)
          .survfit <- summary(survfit(Surv(time,event) ~ grps, data=.data))
          
          .obs <- sapply(1:n.groups, function(x) {
            .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
            .survfit$surv[.indic]
          })
          
          all_pred <- c(all_pred, .est)
          all_obs <- c(all_obs, .obs)
        }
        
        col <- if(hasArg(col)) list(...)$col else 1
        lwd <- if(hasArg(lwd)) list(...)$lwd else 2
        xlim <- if(hasArg(xlim)) list(...)$xlim else c(0,1)
        ylim <- if(hasArg(ylim)) list(...)$ylim else c(0,1)
        xlab <- if(hasArg(xlab)) list(...)$xlab else "Predicted survival"
        ylab <- if(hasArg(ylab)) list(...)$ylab else "Observed survival"
        main <- if(hasArg(main)) list(...)$main else ""
        
        plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
        fit <- loess(all_obs ~ all_pred)
        ord <- order(all_pred)
        lines(all_pred[ord], predict(fit)[ord], col=col, lwd=lwd)
        abline(c(0,1), lty=2)
        
      } else {
      cols <- if (hasArg(col)) rep(list(...)$col, length.out = x$m) else rainbow(x$m)
      ltys <- if (hasArg(lty)) rep(list(...)$lty, length.out = x$m) else rep(1, x$m)
      for (i in 1:x$m) {
        if(is.null(pro.time)){ pro.time <- median(x$data[[i]][,x$outcome[[1]]]) }
        time <- x$data[[i]][,x$outcome[[1]]]
        event <- x$data[[i]][,x$outcome[[2]]]
        .lp = x$calibration[[i]]$lp
        hazard = x$calibration[[i]]$H0.multi
        T.multi = x$calibration[[i]]$time
        
        .pred <- exp(matrix(exp(.lp)) %*% t(as.matrix(-1*hazard)))
        pos = max(which(T.multi <= pro.time))
        .pred = .pred[,pos]
        
        .grps <- as.numeric(cut(.pred,
                                breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                                labels = 1:n.groups))
        
        .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )
        
        .data <- data.frame(time = time, event = event, grps = .grps)
        
        .survfit <- summary(survfit(Surv(time, event) ~ grps, data=.data))
        
        .obs <- sapply(1:n.groups, FUN = function(x) {
          .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
          .survfit$surv[ .indic ] } )
        
        .lower <- sapply(1:n.groups, FUN = function(x) {
          .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
          .survfit$lower[ .indic ] } )
        
        .upper <- sapply(1:n.groups, FUN = function(x) {
          .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
          .survfit$upper[ .indic ] } )
        
        if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
        if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
        if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
        if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
        if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
        if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
        if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
        if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
        if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
        
        if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
        if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}
        
        if(hasArg(ylab)==FALSE) {ylab <- "Observed survival"} else {ylab <- list(...)$ylab}
        if(hasArg(xlab)==FALSE) {xlab <- "Predicted survival"} else {xlab <- list(...)$xlab}
        if(hasArg(main)==FALSE) {main <- ""} else {main <- list(...)$main}
        
        if (i == 1) {
          plot(.est, .obs, type = type, col = cols[i], lty = ltys[i],
               cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
               lwd = lwd, pch = pch, ylim = ylim, xlim = xlim,
               ylab = ylab, xlab = xlab, main = main)
          abline(c(0,1), lty = 2)
        } else {
          points(.est, .obs, type = type, col = cols[i], lty = ltys[i],
                 cex = cex, lwd = lwd, pch = pch)
        }
        
        abline(c(0,1), lty=2)
        
        segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
      }
      }
    } else {
      if(is.null(pro.time)){ pro.time <- median(x$data[,x$outcome[[1]]]) }
      
      if (any(is.na(x$data))){x$data <- na.omit(x$data) }
      
      time <- x$data[,x$outcome[[1]]]
      event <- x$data[,x$outcome[[2]]]
      .lp = x$calibration$lp
      hazard = x$calibration$H0.multi
      T.multi = x$calibration$time
      
      .pred <- exp(matrix(exp(.lp)) %*% t(as.matrix(-1*hazard)))
      pos = max(which(T.multi <= pro.time))
      .pred = .pred[,pos]
      
      .grps <- as.numeric(cut(.pred,
                              breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
      
      .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )
      
      .data <- data.frame(time = time, event = event, grps = .grps)
      
      .survfit <- summary(survfit(Surv(time, event) ~ grps, data=.data))
      
      .obs <- sapply(1:n.groups, FUN = function(x) {
        .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
        .survfit$surv[ .indic ] } )
      
      .lower <- sapply(1:n.groups, FUN = function(x) {
        .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
        .survfit$lower[ .indic ] } )
      
      .upper <- sapply(1:n.groups, FUN = function(x) {
        .indic <- sum(as.numeric(.survfit$strata)==x & .survfit$time<=pro.time)
        .survfit$upper[ .indic ] } )
      
      if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
      if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
      if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
      if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
      if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
      if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
      if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
      if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
      if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
      
      if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
      if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}
      
      if(hasArg(ylab)==FALSE) {ylab <- "Observed survival"} else {ylab <- list(...)$ylab}
      if(hasArg(xlab)==FALSE) {xlab <- "Predicted survival"} else {xlab <- list(...)$xlab}
      if(hasArg(main)==FALSE) {main <- ""} else {main <- list(...)$main}
      
      plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
           type = type, col = col, lty = lty, lwd = lwd, main=main,
           pch = pch, ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)
      
      abline(c(0,1), lty=2)
      
      segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
    }
    
  }
}