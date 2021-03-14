#' Print method for plsRglm models
#' 
#' This function provides a predict method for the class \code{"plsRglmmodel"}
#' 
#' 
#' @param object An object of the class \code{"plsRmodel"}.
#' @param newdata An optional data frame in which to look for variables with
#' which to predict. If omitted, the fitted values are used.
#' @param comps A value with a single value of component to use for prediction.
#' @param type Type of predicted value. Available choices are the glms ones
#' ("\code{link}", "\code{response}", "\code{terms}"), the polr ones
#' ("\code{class}", "\code{probs}") or the scores ("\code{scores}").
#' @param se.fit If TRUE, pointwise standard errors are produced for the
#' predictions using the Cox model.
#' @param weights Vector of case weights. If \code{weights} is a vector of
#' integers, then the estimated coefficients are equivalent to estimating the
#' model from data with the individual \code{cases} replicated as many times as
#' indicated by \code{weights}.
#' @param dispersion the dispersion of the GLM fit to be assumed in computing
#' the standard errors. If omitted, that returned by summary applied to the
#' object is used.
#' @param methodNA Selects the way of predicting the response or the scores of
#' the new data. For complete rows, without any missing value, there are two
#' different ways of computing the prediction. As a consequence, for mixed
#' datasets, with complete and incomplete rows, there are two ways of computing
#' prediction : either predicts any row as if there were missing values in it
#' (\code{missingdata}) or selects the prediction method accordingly to the
#' completeness of the row (\code{adaptative}).
#' @param verbose should info messages be displayed ?
#' @param \dots Arguments to be passed on to \code{stats::glm} and
#' \code{plsRglm::plsRglm}.
#' @return When type is "\code{response}", a matrix of predicted response
#' values is returned.\cr When type is "\code{scores}", a score matrix is
#' returned.
#' @author Frédéric Bertrand\cr
#' \email{frederic.bertrand@@math.unistra.fr}\cr
#' \url{https://fbertran.github.io/homepage/}
#' @seealso See Also \code{\link[stats]{predict.glm}}
#' @references Nicolas Meyer, Myriam Maumy-Bertrand et
#' Frédéric Bertrand (2010). Comparing the linear and the
#' logistic PLS regression with qualitative predictors: application to
#' allelotyping data. \emph{Journal de la Societe Francaise de Statistique},
#' 151(2), pages 1-18.
#' \url{http://publications-sfds.math.cnrs.fr/index.php/J-SFdS/article/view/47}
#' @keywords methods predict
#' @examples
#' 
#' data(pine)
#' Xpine<-pine[,1:10]
#' ypine<-pine[,11]
#' data(pine_sup)
#' Xpine_sup<-pine_sup[,1:10]
#' Xpine_supNA<-Xpine_sup
#' Xpine_supNA[1,1]<-NA
#' 
#' modpls=plsRglm(dataY=ypine,dataX=Xpine,nt=6,modele="pls-glm-family",family="gaussian",
#' verbose=FALSE)
#' modplsform=plsRglm(x11~.,data=pine,nt=6,modele="pls-glm-family",family="gaussian", verbose=FALSE)
#' modpls2=plsRglm(dataY=ypine,dataX=Xpine,nt=6,modele="pls-glm-family",
#' dataPredictY=Xpine_sup,family="gaussian", verbose=FALSE)
#' modpls2NA=plsRglm(dataY=ypine,dataX=Xpine,nt=6,modele="pls-glm-family",
#' dataPredictY=Xpine_supNA,family="gaussian", verbose=FALSE)
#' 
#' #Identical to predict(modpls,type="link") or modpls$Std.ValsPredictY
#' cbind(modpls$Std.ValsPredictY,modplsform$Std.ValsPredictY,
#' predict(modpls),predict(modplsform))
#' 
#' #Identical to predict(modpls,type="response") or modpls$ValsPredictY
#' cbind(modpls$ValsPredictY,modplsform$ValsPredictY,
#' predict(modpls,type="response"),predict(modplsform,type="response"))
#' 
#' #Identical to modpls$ttPredictY
#' predict(modpls,type="scores")
#' predict(modplsform,type="scores")
#' 
#' \donttest{
#' #Identical to modpls2$ValsPredictY
#' cbind(predict(modpls,newdata=Xpine_sup,type="response"),
#' predict(modplsform,newdata=Xpine_sup,type="response"))
#' 
#' #Select the number of components to use to derive the prediction
#' predict(modpls,newdata=Xpine_sup,type="response",comps=1)    
#' predict(modpls,newdata=Xpine_sup,type="response",comps=3)    
#' predict(modpls,newdata=Xpine_sup,type="response",comps=6)    
#' try(predict(modpls,newdata=Xpine_sup,type="response",comps=8))
#' 
#' #Identical to modpls2$ttValsPredictY
#' predict(modpls,newdata=Xpine_sup,type="scores")    
#' 
#' #Select the number of components in the scores matrix
#' predict(modpls,newdata=Xpine_sup,type="scores",comps=1)    
#' predict(modpls,newdata=Xpine_sup,type="scores",comps=3)    
#' predict(modpls,newdata=Xpine_sup,type="scores",comps=6)    
#' try(predict(modpls,newdata=Xpine_sup,type="scores",comps=8))
#' 
#' #Identical to modpls2NA$ValsPredictY
#' predict(modpls,newdata=Xpine_supNA,type="response",methodNA="missingdata")    
#' 
#' cbind(predict(modpls,newdata=Xpine_supNA,type="response"),
#' predict(modplsform,newdata=Xpine_supNA,type="response"))
#' 
#' predict(modpls,newdata=Xpine_supNA,type="response",comps=1)    
#' predict(modpls,newdata=Xpine_supNA,type="response",comps=3)    
#' predict(modpls,newdata=Xpine_supNA,type="response",comps=6)    
#' try(predict(modpls,newdata=Xpine_supNA,type="response",comps=8))
#' 
#' #Identical to modpls2NA$ttPredictY
#' predict(modpls,newdata=Xpine_supNA,type="scores",methodNA="missingdata")
#' predict(modplsform,newdata=Xpine_supNA,type="scores",methodNA="missingdata")
#' 
#' predict(modpls,newdata=Xpine_supNA,type="scores")    
#' predict(modplsform,newdata=Xpine_supNA,type="scores")    
#' predict(modpls,newdata=Xpine_supNA,type="scores",comps=1)    
#' predict(modpls,newdata=Xpine_supNA,type="scores",comps=3)    
#' predict(modpls,newdata=Xpine_supNA,type="scores",comps=6)    
#' try(predict(modpls,newdata=Xpine_supNA,type="scores",comps=8))
#' }    
#' 
#' @export
predict.plsRglmmodel <- function(object,newdata,comps=object$computed_nt,type=c("link", "response", "terms", "scores", "class", "probs"),se.fit=FALSE,weights, dispersion = NULL,methodNA="adaptative",verbose=TRUE,...)
{
  if (!inherits(object, "plsRglmmodel")) 
    stop("Primary argument much be a plsRglmmodel object")
  if(missing(type)){type="link"}
  type <- match.arg(type)
  if (!(type %in% c("link", "response", "terms", "scores", "class", "probs")))
    stop("Invalid type specification")
  if (comps>object$computed_nt){stop("Cannot predict using more components than extracted.")}
  if (missing(newdata) || is.null(newdata)) {
    nrtt <- nrow(object$tt)
    if (type=="link"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
                      if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "link",se.fit=se.fit, dispersion = dispersion,...))}
    }
    if (type=="response"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
                      if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit, dispersion = dispersion,...))}                      
    }
    if (type=="terms"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
                       if("glm" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit, dispersion = dispersion,...))}                       
    }  
    if (type=="class"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
                       if("polr" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "class",...))}                       
    }  
    if (type=="probs"){ttpredY<-data.frame(cbind(object$tt[,1:comps],matrix(0,nrow=nrtt,ncol=object$computed_nt-comps)));colnames(ttpredY) <- NULL
                       if("polr" %in% class(object$FinalModel)){ttpred<-data.frame(tt=ttpredY);return(predict(object$FinalModel,newdata=ttpred,type = "probs",...))}                       
    }  
    if (type=="scores"){return(object$tt[,1:comps])}   
  } else {
    nrnd <- nrow(newdata)
    if(any(apply(is.na(newdata),MARGIN=1,"all"))){return(vector("list",0)); cat("One of the rows of newdata is completely filled with missing data\n"); stop()}
    if(any(is.na(newdata))) {na.miss.newdata <- TRUE} else {na.miss.newdata <- FALSE}
    if(!is.null(object$call$formula)){
      mf <- match.call(expand.dots = FALSE)
      m <- match(c("subset", "weights"), names(mf), 0L)
      mf <- mf[c(1L, m)]
      mf$data <- newdata
      mf$formula <- object$call$formula[-2]
      mf$drop.unused.levels <- TRUE
      mf$na.action <- na.pass    
      mf[[1L]] <- as.name("model.frame")
      mf <- eval(mf, parent.frame())
      mt <- attr(mf, "terms")
      #    attr(mt,"intercept")<-0L    
      newdata.frame <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts)[,-1]
      else matrix(, nrnd, 0L)
      weights <- as.vector(model.weights(mf))
    } else {newdata.frame <- newdata}
    newdata.scaled <- sweep(sweep(newdata.frame, 2, attr(object$ExpliX,"scaled:center")), 2 ,attr(object$ExpliX,"scaled:scale"), "/")
  newdataNA <- !is.na(newdata)
  newdata.scaledwotNA <- as.matrix(newdata.scaled)
  newdata.scaledwotNA[!newdataNA] <- 0
  ttpredY <- NULL
  if (methodNA=="adaptative") {
    for(ii in 1:nrnd){
      if (all(newdataNA[ii,])){
        ttpredY <- rbind(ttpredY, c(newdata.scaledwotNA[ii,]%*%object$wwetoile[,1:comps],rep(0,object$computed_nt-comps))) 
      }
      else {
        if(verbose){cat("Missing value in row ",ii,".\n")}
        ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
      }}}
  if(methodNA=="missingdata") {
    if(verbose){cat("Prediction as if missing values in every row.\n")}
    for (ii in 1:nrnd) {  
      ttpredY <- rbind(ttpredY, c(t(solve(t(object$pp[newdataNA[ii,],,drop=FALSE])%*%object$pp[newdataNA[ii,],,drop=FALSE])%*%t(object$pp[newdataNA[ii,],,drop=FALSE])%*%(newdata.scaledwotNA[ii,])[newdataNA[ii,]])[1:comps],rep(0,object$computed_nt-comps)))
    }
  }
colnames(ttpredY) <- NULL
ttpred<-data.frame(tt=ttpredY)
  if("glm" %in% class(object$FinalModel)){
  if (type=="link"){return(predict(object$FinalModel,newdata=ttpred,type = "link",se.fit=se.fit,dispersion = dispersion,...))                  
  }
  if (type=="response"){return(predict(object$FinalModel,newdata=ttpred,type = "response",se.fit=se.fit,dispersion = dispersion,...))                    
  }
  if (type=="terms"){return(predict(object$FinalModel,newdata=ttpred,type = "terms",se.fit=se.fit,dispersion = dispersion,...))
  }
  }
  if("polr" %in% class(object$FinalModel)){
  if (type=="class"){return(predict(object$FinalModel,newdata=ttpred,type = "class",...))}  
  if (type=="probs"){return(predict(object$FinalModel,newdata=ttpred,type = "probs",...))}  
  }
  if (type=="scores"){colnames(ttpredY) <- paste("Comp_",1:object$computed_nt,sep="");
                      return(ttpredY[,1:comps,drop=FALSE])}
}
}

