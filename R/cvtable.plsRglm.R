cvtable.plsRglm <- function(x,verbose=TRUE,...)
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
    mincvMC<-function(lll){safe_which_min(lll[-1,4])} 
    mincvMCobs<-sapply(x,mincvMC)
    rescvMC<-count_levels(mincvMCobs, start_level = 1L) 
    if(verbose){print(rescvMC)}
  }
  
  if("Q2Chisq_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV Q2Chi2 criterion:")}
    mincvQ2Chisq<-function(lll){threshold_crossing(lll[-1,5+2*MClassed], lll[-1,4+2*MClassed])}  
    mincvQ2Chisqobs<-sapply(x,mincvQ2Chisq)
    rescvQ2Chisq<-count_levels(mincvQ2Chisqobs, start_level = 0L)
    if(verbose){print(rescvQ2Chisq)}    
    }
  
  if("PREChi2_Pearson_Y" %in% colnames(x[[1]])){
    if(verbose){cat("\nCV PreChi2 criterion:")}
    mincvPreChi2<-function(lll){safe_which_min(lll[-1,6+2*MClassed])} 
    mincvPreChi2obs<-sapply(x,mincvPreChi2)     
    rescvPreChi2<-count_levels(mincvPreChi2obs, start_level = 1L)
    if(verbose){print(rescvPreChi2)}
  }
  
  if(MClassed){
    res=list(CVMC=rescvMC,CVQ2Chi2=rescvQ2Chisq,CVPreChi2=rescvPreChi2)
    class(res) <- "table.summary.cv.plsRglmmodel"
    invisible(res)    
  } else {
    res=list(CVQ2Chi2=rescvQ2Chisq,CVPreChi2=rescvPreChi2)  
    class(res) <- "table.summary.cv.plsRglmmodel"
    invisible(res)    
  }  
}
