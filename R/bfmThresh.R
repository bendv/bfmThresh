#' @title BFM With Magnitude Threshold
#' @description bfm with automatic magn thresholding and re-monitoring
#' 
#' @param x Numeric. Time series vector
#' @param dates Time (date) index of \code{x}
#' @param start Numeric. Start of monitoring period. See \code{\link{bfastmonitor}}
#' @param monend Numeric. Optional: end of the monitoring period. See \code{\link{bfmSpatial}}
#' @param magnThresh Numeric. Threshold for breakpoints. All breakpoints with magnitude < \code{magnThresh} are ignored, all data \code{adjust} years preceding this breakpoint are deleted, and monitoring continues with a revised history period.
#' @param magnPeriod Numeric. Period of time (in years) from which to compute magnitude (as median of residuals)
#' @param adjust Numeric. Delete all observations within this many years preceding a breakpoint, if that breakpoint has a magnitude less than \code{magnThresh}
#' @param verbose Logical. Display error messages (if any)?
#' @param formula See \code{\link{bfastmonitor}}
#' @param order See \code{\link{bfastmonitor}}
#' @param returnBFM Logical. Return a regular \code{bfastmonitor} object?
#' @param ... Arguments to pass to \code{\link{bfastmonitor}}
#' 
#' @return Numeric vector with output: break date, magnitude, length of history period (in years), R-squared of fitted model, Adjusted R-squared of fitted model, Intercept of fitted model, cosine coefficient(s) (if applicable), sine coefficient(s) (if applicable), error flag (1 means an error was encountered in the \code{try()} statement). If \code{returnBFM} is \code{FALSE}, a list with the above values as first item and a regular BFM object as second item is returned.
#' 
#' @details
#' The magnitude is computed as the median of the residuals within the period \code{magnPeriod} years after the break date. If this magnitude is greater than \code{magnThresh}, all observations preceding the break date by \code{adjust} years are deleted (including the observation at the break date), the start of the monitoring period is reset to the break date, the model is re-fit and monitoring continues until the end of the time series is reached, or a breakpoint with magnitude > \code{magnThresh} is found.
#' 
#' The default value for \code{magnPeriod} is 1, assuming that no follow-up land use or other dynamics take place within a 1-year period.
#' 
#' The amplitude of the fitted 1st order model can be computed by taking \code{sqrt{harmoncos^2 + harmonsin^2}}
#' 
#' @author Ben DeVries
#' @import bfastSpatial
#' 
#' @examples
#' library(bfastSpatial)
#' data(tura)
#' s <- getSceneinfo(names(tura))
#' xy <- c(820692.4, 829564.4)
#' z <- zoo(tura[cellFromXY(tura, xy)][1, ], s$date)
#' z <- window(z, start = "1999-01-01")
#' bfm <- bfmThresh(z, dates = time(z), start = c(2005, 1), plot = TRUE)
#' print(bfm)
#' bfm2 <- bfmThresh(z, dates = time(z), start = c(2005, 1), magnThresh = -500, plot = TRUE)
#' print(bfm2)


bfmThresh <- function(x, dates, start, monend = NULL, magnThresh = 0, magnPeriod = 1, adjust = 1, verbose = FALSE, formula = response ~ harmon, order = 1, returnBFM = FALSE, ...) {
  bts <- bfastts(x, dates, type = 'irregular')
  
  # trim time series by monend
  if(!is.null(monend))
    bts <- window(bts, end = monend)
  
  # magnitude monitor
  iterate <- TRUE
  err <- NA
  
  # determine length of coefficient vector
  # = intercept [+ trend] [+ harmoncos*order] [+ harmonsin*order]
  coef_len <- 1 # intercept
  modterms <- attr(terms(formula), "term.labels")
  if("trend" %in% modterms)
    coef_len <- coef_len + 1
  if("harmon" %in% modterms)
    coef_len <- coef_len + (order * 2) # sin and cos terms
  
  while(iterate) {
    bfm <- try(bfastmonitor(bts, start = start, formula = formula, order = order, ...), silent = !verbose)
    if(class(bfm) == 'try-error') {
      bkp <- magn <- hist <- rsq <- adjrsq <- NA
      coefs <- rep(NA, coef_len)
      err <- 1
      break
    } else {
      bkp <- bfm$breakpoint
      hist <- bfm$history[2] - bfm$history[1]
      rsq <- summary(bfm$model)$r.squared
      adjrsq <- summary(bfm$model)$adj.r.squared
      coefs <- coef(bfm$model)
    }
    
    # breakpoint magnitude
    if(!is.na(bkp)) {
      postpp <- subset(bfm$tspp, time >= bkp & time <= (bkp + magnPeriod))
      magn <- median(postpp$response - postpp$prediction)
    } else {
      magn <- NA
    }
    
    if(!is.na(magn) & magn >= magnThresh) {
      bts[(time(bts) <= bkp) & (time(bts) > bkp - adjust)] <- NA
      start <- c(floor(bkp), 1)
    } else {
      iterate <- FALSE
    }
  }
  
  res <- c(bkp, magn, hist, rsq, adjrsq, coefs, err)
  names(res) <- c("breakpoint", "magnitude", "history", "r.squared", "adj.r.squared", names(coefs), "error")
  
  if(returnBFM)
    res <- list(res, bfm)
  
  return(res)
  
}