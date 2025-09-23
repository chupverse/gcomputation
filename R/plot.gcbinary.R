plot.gcbinary <- function (x, method="calibration", n.groups=5, smooth=FALSE, ...) {
  if (!(method %in% c("calibration","proportion"))) {stop("Method needs to be calibration or proportion")}
  if (!is.null(x$newdata)) {stop("Plots do not work on a transposed object, use the original object")}
  if (method == "proportion") {
    data = x$data
    outcome = x$outcome
    group = x$group
    
    datag0 = data[which(data[,group] == 0),]
    datag1 = data[which(data[,group] == 1),]
    
    #### This part unlike survival model calibration same as obs
    predict.p0 = mean(x$calibration$predict[which(data[,group] == 0)])
    predict.p1 = mean(x$calibration$predict[which(data[,group] == 1)])
    obs.p0 = nrow(datag0[which(datag0[,outcome] == 1),]) / length(datag0[,outcome])
    obs.p1 = nrow(datag1[which(datag1[,outcome] == 1),]) / length(datag1[,outcome])
    
    
    if (hasArg(labels) == FALSE) {
      x_labels <- c("Group 0:\nObserved", "Group 0:\nPredicted", "Group 1:\nObserved", "Group 1:\nPredicted")
    } else {
      x_labels <- list(...)$labels
    }
    
    plot(x=c(1,1.2,2,2.2), y=c(obs.p0,predict.p0,obs.p1,predict.p1),xlab="",ylab='',main=NULL, xaxt="n")
    axis(1, at = c(1, 1.2, 2, 2.2), labels = x_labels)
    
  }
  

  
  if (method == "calibration") {
    
    if (!is.null(x$m)) {
      
      if (smooth) {
        all_pred <- unlist(lapply(1:x$m, function(i) x$calibration[[i]]$predict))
        all_outcome <- unlist(lapply(1:x$m, function(i) x$data[[i]][, x$outcome]))
        
        col <- if(hasArg(col)) list(...)$col else 1
        lwd <- if(hasArg(lwd)) list(...)$lwd else 2
        xlim <- if(hasArg(xlim)) list(...)$xlim else c(0,1)
        ylim <- if(hasArg(ylim)) list(...)$ylim else c(0,1)
        xlab <- if(hasArg(xlab)) list(...)$xlab else "Predicted proportions"
        ylab <- if(hasArg(ylab)) list(...)$ylab else "Observed proportions"
        main <- if(hasArg(main)) list(...)$main else ""
        
        plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
        
        fit <- loess(all_outcome ~ all_pred)
        ord <- order(all_pred)
        lines(all_pred[ord], predict(fit)[ord], col=col, lwd=lwd)
        abline(c(0,1), lty=2)
        
      } else {
      
      cols <- if (hasArg(col)) rep(list(...)$col, length.out = x$m) else rainbow(x$m)
      ltys <- if (hasArg(lty)) rep(list(...)$lty, length.out = x$m) else rep(1, x$m)
      for (i in 1:x$m) {
        .event = x$data[[i]][,x$outcome]
        .pred = x$calibration[[i]]$predict
        data = x$data[[i]]
        outcome = x$data[[i]][,x$outcome]
        
        if (length(unique(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))))) != length(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))))) {
          stop("n.groups too high")
        }
        
        .grps <- as.numeric(cut(.pred,
                                breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                                labels = 1:n.groups))
        
        .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )
        
        .data <- data.frame(outcome = outcome, grps = .grps)
        
        .mod <- glm(outcome ~ grps, data=.data)
        .fit <- summary(glm(outcome ~ grps, data=.data))
        
        .obs <- sapply(1:n.groups, FUN = function(x) {
          .indic <- .mod$data$grps == x
          mean(.mod$data$outcome[.indic])} )
        
        .lower <- sapply(1:n.groups, function(x) {
          .indic <- .mod$data$grps == x
          p <- mean(.mod$data$outcome[.indic])
          n <- sum(.indic)
          p - 1.96 * sqrt(p * (1 - p) / n)
        })
        
        .upper <- sapply(1:n.groups, function(x) {
          .indic <- .mod$data$grps == x
          p <- mean(.mod$data$outcome[.indic])
          n <- sum(.indic)
          p + 1.96 * sqrt(p * (1 - p) / n)
        })
        
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
        
        if(hasArg(ylab)==FALSE) {ylab <- "Observed proportions"} else {ylab <- list(...)$ylab}
        if(hasArg(xlab)==FALSE) {xlab <- "Predicted proportions"} else {xlab <- list(...)$xlab}
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
      
      if (any(is.na(x$data))){x$data <- na.omit(x$data) }
      
      .event = x$data[,x$outcome]
      .pred = x$calibration$predict
      data = x$data
      outcome = x$data[,x$outcome]
      
      if (length(unique(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))))) != length(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))))) {
        stop("n.groups too high")
      }
      
      .grps <- as.numeric(cut(.pred,
                              breaks = c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
      
      .est <- sapply(1:n.groups, FUN = function(x) { mean(.pred[.grps==x]) } )
      
      .data <- data.frame(outcome = outcome, grps = .grps)
      
      .mod <- glm(outcome ~ grps, data=.data)
      .fit <- summary(glm(outcome ~ grps, data=.data))
      
      .obs <- sapply(1:n.groups, FUN = function(x) {
        .indic <- .mod$data$grps == x
        mean(.mod$data$outcome[.indic])} )
      
      .lower <- sapply(1:n.groups, function(x) {
        .indic <- .mod$data$grps == x
        p <- mean(.mod$data$outcome[.indic])
        n <- sum(.indic)
        p - 1.96 * sqrt(p * (1 - p) / n)
      })
      
      .upper <- sapply(1:n.groups, function(x) {
        .indic <- .mod$data$grps == x
        p <- mean(.mod$data$outcome[.indic])
        n <- sum(.indic)
        p + 1.96 * sqrt(p * (1 - p) / n)
      })
      
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
      
      if(hasArg(ylab)==FALSE) {ylab <- "Observed proportions"} else {ylab <- list(...)$ylab}
      if(hasArg(xlab)==FALSE) {xlab <- "Predicted proportions"} else {xlab <- list(...)$xlab}
      if(hasArg(main)==FALSE) {main <- ""} else {main <- list(...)$main}
      
      plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
           type = type, col = col, lty = lty, lwd = lwd, main=main,
           pch = pch, ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)
      
      abline(c(0,1), lty=2)
      
      segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
      
    }
    
    
  }
}
