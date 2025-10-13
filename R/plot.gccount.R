plot.gccount <- function (x, method="calibration", n.groups=5, smooth=FALSE, ...) {
  if (!(method %in% c("calibration","proportion"))) {stop("Method needs to be calibration or proportion")}
  if (!is.null(x$newdata)) {stop("Plots do not work on a transposed object, use the original object")}
  if (method == "proportion") {
    if (!is.null(x$m)) {stop("The \"method=proportion\" is not available when \"boot.mi=TRUE\"")}
    
    data = x$data
    outcome = x$outcome
    group = x$group
    
    datag0 = data[which(data[,group] == 0),]
    datag1 = data[which(data[,group] == 1),]
    
    #### This part unlike survival model calibration same as obs
    predict.c0 = mean(x$calibration$predict[which(data[,group] == 0)])
    predict.c1 = mean(x$calibration$predict[which(data[,group] == 1)])
    obs.c0 = mean(datag0[[outcome]]) 
    obs.c1 = mean(datag1[[outcome]])
    
    
    if (hasArg(labels) == FALSE) {
      x_labels <- c("Group 0:\nObserved", "Group 0:\nPredicted", "Group 1:\nObserved", "Group 1:\nPredicted")
    } else {
      x_labels <- list(...)$labels
    }
    
    plot(x=c(1,1.2,2,2.2), y=c(obs.c0,predict.c0,obs.c1,predict.c1),xlab="",ylab='',main=NULL, xaxt="n")
    axis(1, at = c(1, 1.2, 2, 2.2), labels = x_labels)
    
  }
  
  
  
  if (method == "calibration") {
    
    if (!is.null(x$m)) {
      cols <- if (hasArg(col)) rep(list(...)$col, length.out = x$m) else rainbow(x$m)
      ltys <- if (hasArg(lty)) rep(list(...)$lty, length.out = x$m) else rep(1, x$m)
      n.grouperror <- 0
      firstplot <- TRUE
      all_est = c()
      all_obs = c()
      all_lower = c()
      all_upper = c()
      for (i in 1:x$m) {
        .pred = x$calibration[[i]]$predict
        data = x$data[[i]]
        outcome = x$data[[i]][,x$outcome]
        
        if (length(unique(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))))) != length(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))))) {
          n.grouperror <- n.grouperror + 1
          next
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
          l <- mean(.mod$data$outcome[.indic])
          n <- sum(.indic)
          l - 1.96 * sqrt(l/n)
        })
        
        .upper <- sapply(1:n.groups, function(x) {
          .indic <- .mod$data$grps == x
          l <- mean(.mod$data$outcome[.indic])
          n <- sum(.indic)
          l + 1.96 * sqrt(l/n)
        })
        
        if (any(sapply(.lower, function(x) length(x) == 0))) {
          n.grouperror <- n.grouperror + 1
          next
        }
        
        if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
        if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
        if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
        if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
        if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
        if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
        if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
        if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
        if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
        
        if(hasArg(ylim)==FALSE) {ylim <- c(0,max(outcome))} else {ylim <- list(...)$ylim}
        if(hasArg(xlim)==FALSE) {xlim  <- c(0,max(outcome))} else {xlim <- list(...)$xlim}
        
        if(hasArg(ylab)==FALSE) {ylab <- "Observed events number"} else {ylab <- list(...)$ylab}
        if(hasArg(xlab)==FALSE) {xlab <- "Predicted events number"} else {xlab <- list(...)$xlab}
        if(hasArg(main)==FALSE) {main <- ""} else {main <- list(...)$main}
        
        
        if (smooth == TRUE) {
          all_est = c(all_est, .est)
          all_obs = c(all_obs, .obs)
          all_lower = c(all_lower, .lower)
          all_upper = c(all_upper, .upper)
        } else {
          if (firstplot == TRUE) {
            plot(.est, .obs, type = type, col = cols[i], lty = ltys[i],
                 cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
                 lwd = lwd, pch = pch, ylim = ylim, xlim = xlim,
                 ylab = ylab, xlab = xlab, main = main)
            firstplot = FALSE
          } else {
            points(.est, .obs, type = type, col = cols[i], lty = ltys[i],
                   cex = cex, lwd = lwd, pch = pch)
          }
          abline(c(0,1), lty=2)
          segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
        }
      }
      if (n.grouperror > 0) {
        if (n.grouperror == 1) {
          warning(paste0("There was ",n.grouperror," model where \"n.groups\" was too high (See ?plot.gcbinary)"))
        } else {
          warning(paste0("There were ",n.grouperror," models where \"n.groups\" was too high (See ?plot.gcbinary)"))
        }
      }
      if (smooth == TRUE) {
        plot(0, 0, type="n", xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main)
        
        loesssmooth <- loess(all_obs ~ all_est)
        loesslower <- loess(all_lower ~ all_est)
        loessupper <- loess(all_upper ~ all_est)
        
        all_est_grid <- seq(min(all_est), max(all_est), length.out=200)
        lines(all_est_grid, predict(loesssmooth, newdata=data.frame(all_est=all_est_grid)), col=col, lwd=lwd, lty=lty)
        lines(all_est_grid, predict(loesslower, newdata=data.frame(all_est=all_est_grid)), col=col, lty=2)
        lines(all_est_grid, predict(loessupper, newdata=data.frame(all_est=all_est_grid)), col=col, lty=2)
        
        abline(c(0,1), lty=2)
      }
      
    } else {
      
      if (any(is.na(x$data))){x$data <- na.omit(x$data) }
      .pred = x$calibration$predict
      data = x$data
      outcome = x$data[,x$outcome]
      
      print(.pred)
      print(unique(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups)))))
      print(c(-Inf, quantile(.pred, seq(1/n.groups, 1, 1/n.groups))))
      
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
        l <- mean(.mod$data$outcome[.indic])
        n <- sum(.indic)
        l - 1.96 * sqrt(l/n)
      })
      
      .upper <- sapply(1:n.groups, function(x) {
        .indic <- .mod$data$grps == x
        l <- mean(.mod$data$outcome[.indic])
        n <- sum(.indic)
        l + 1.96 * sqrt(l/n)
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
      
      if(hasArg(ylim)==FALSE) {ylim <- c(0,max(outcome))} else {ylim <- list(...)$ylim}
      if(hasArg(xlim)==FALSE) {xlim  <- c(0,max(outcome))} else {xlim <- list(...)$xlim}
      
      if(hasArg(ylab)==FALSE) {ylab <- "Observed events number"} else {ylab <- list(...)$ylab}
      if(hasArg(xlab)==FALSE) {xlab <- "Predicted events number"} else {xlab <- list(...)$xlab}
      if(hasArg(main)==FALSE) {main <- ""} else {main <- list(...)$main}
      
      plot(.est, .obs, cex = cex, cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
           type = type, col = col, lty = lty, lwd = lwd, main=main,
           pch = pch, ylim = ylim, xlim = xlim, ylab=ylab, xlab=xlab)
      
      abline(c(0,1), lty=2)
      
      segments(x0 = .est, y0 = .lower, x1 = .est, y1 = .upper, col = col, lwd = lwd)
      
    }
    
    
  }
}
