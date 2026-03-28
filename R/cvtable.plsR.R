cvtable.plsR <- function(x,verbose=TRUE,...)
{
  safe_which_min <- function(values) {
    values <- as.numeric(values)
    finite <- is.finite(values)
    if (!any(finite)) {
      return(NA_integer_)
    }
    which.min(replace(values, !finite, Inf))
  }

  threshold_crossing <- function(metric, limit) {
    metric <- as.numeric(metric)
    limit <- as.numeric(limit)
    valid_idx <- which(is.finite(metric) & is.finite(limit))
    if (!length(valid_idx)) {
      return(NA_integer_)
    }
    failed_idx <- valid_idx[!(metric[valid_idx] > limit[valid_idx])]
    if (!length(failed_idx)) {
      return(valid_idx[length(valid_idx)])
    }
    failed_idx[1] - 1L
  }

  count_levels <- function(obs, start_level = 0L) {
    obs <- as.integer(obs)
    valid_obs <- obs[!is.na(obs)]
    max_level <- if (length(valid_obs)) max(valid_obs) else start_level
    table(factor(valid_obs, levels = seq.int(start_level, max_level)))
  }

  MClassed=FALSE
  if("CV_MissClassed" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV MissClassed criterion:")}
    MClassed=TRUE
    mincvMC<-function(lll){safe_which_min(lll[-1,3])} 
    mincvMCobs<-sapply(x,mincvMC)
    rescvMC<-count_levels(mincvMCobs, start_level = 1L) 
    if(verbose){print(rescvMC)}
  }
  
  if("Q2_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV Q2 criterion:")}
    mincvQ2<-function(lll){threshold_crossing(lll[-1,4+2*MClassed], lll[-1,3+2*MClassed])}  
    mincvQ2obs<-sapply(x,mincvQ2)
    rescvQ2<-count_levels(mincvQ2obs, start_level = 0L)   
    if(verbose){print(rescvQ2)}    
  }
  
  if("PRESS_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV Press criterion:")}
    mincvPress<-function(lll){safe_which_min(lll[-1,5+2*MClassed])} 
    mincvPressobs<-sapply(x,mincvPress)     
    rescvPress<-count_levels(mincvPressobs, start_level = 1L)
    if(verbose){print(rescvPress)}
  }
  
  if(MClassed){
    res=list(CVMC=rescvMC,CVQ2=rescvQ2,CVPress=rescvPress)
    class(res) <- "table.summary.cv.plsRmodel"
  } else {
    res=list(CVQ2=rescvQ2,CVPress=rescvPress)  
    class(res) <- "table.summary.cv.plsRmodel"
  }

  if (inherits(x, "summary.cv.plsRmultiModel")) {
    class(res) <- c("table.summary.cv.plsRmultiModel", class(res))
  }

  invisible(res)
}
