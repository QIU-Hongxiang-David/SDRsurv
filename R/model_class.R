#' @title S3 class for predictive models with an intercept only
#' @name intercept_model
#' @param value the predicted mean value
#' @return an "`intercept_model`" object, essentially a list with element `value`.
#' @export
intercept_model<-function(value){
    assert_that(is.number(value))
    out<-list(value=value)
    
    class(out)<-"intercept_model"
    out
}

#' @title S3 class for predictive models with an intercept only and influence function information
#' @name intercept_IF_model
#' @param value the predicted mean value
#' @param IF named vector of influence function in one time window evaluated at each observation. The name corresponds to each observation's id. Used to calculate standard error and confidence interval for SDR estimator of a scalar estimand (rather than a function).
#' @return an "`intercept_IF_model`" object, essentially a list with element `value` and `IF`.
#' @export
intercept_IF_model<-function(value,IF){
    assert_that(is.number(value))
    assert_that(is.vector(IF,mode="numeric"))
    assert_that(!is.null(names(IF)))
    out<-list(value=value,IF=IF)
    
    class(out)<-c("intercept_IF_model")
    out
}

#' @export
predict.intercept_model<-function(object,...){
    object$value
}

#' @export
fitted.intercept_IF_model<-predict.intercept_model
#' @export
predict.intercept_IF_model<-predict.intercept_model
#' @export
fitted.intercept_model<-predict.intercept_model


est_IF_SE_CI_conflevel<-function(est,IF,SE,CI,conflevel){
    out<-list(est=est,SE=SE,IF=IF,CI=CI,conf.level=conflevel)
    class(out)<-"est_IF_SE_CI_conflevel"
    out
}

#' @export
print.est_IF_SE_CI_conflevel<-function(est_IF_SE_CI_conflevel.obj,...){
    cat("\n")
    cat("Point estimate: ",est_IF_SE_CI_conflevel.obj$est,"\n")
    cat("Standard error: ",est_IF_SE_CI_conflevel.obj$SE,"\n")
    cat(round(est_IF_SE_CI_conflevel.obj$conf.level*100,2),"% confidence interval:\n")
    print(est_IF_SE_CI_conflevel.obj$CI)
}


#' @title S3 class for the predictive model for (conditional) survival probability with multiple stages
#' @name mult_stage_survfit
#' @description
#' This function is the creator of `mult_stage_survfit` class. A `mult_stage_survfit` object is a collection of fitted models for conditional survival probability within each time window. The associated `predict` and `fitted` methods autocmatically conduct the prediction, which is multiplying all model predictions, potentially after clipping the predictions into an interval in \eqn{[0,1]}.
#' @param covariate.data data frame containing the covariates used to train the models
#' @param formula formula used in the models
#' @param visit.times numeric/integer vector of visit times in ascending order. The first visit time is typically the baseline.
#' @param tvals time t for which P(T > t) (potentially given covariates) is estimated (T is the time to event)
#' @param truncation.index index of the visit time to which left-truncation is applied. The truncation time is `visit.times[truncation.index]`. Covariates available up to (inclusive) `visit.times[truncation.index]` are of interest.
#' @param models list of models. Each model can predict the survival probability beyond the next visit time given survival at the current visit time. Models should be in an ascending order according to the corresponding visit times.
#' @return a "`mult_stage_survfit`" object, essentially a list with input parameters of this function.
#' @export
mult_stage_survfit<-function(covariate.data,formula,visit.times,tval,truncation.index,models){
    K<-length(visit.times)
    
    assert_that(is.list(models))
    assert_that(is_formula(formula))
    
    #check if visit.times are ascending with unique values
    if(is.unsorted(visit.times,strictly=TRUE)){
        stop("visit.times is not sorted in ascending order with unique values")
    }
    
    assert_that(is.number(tval),noNA(tval))
    if(tval<visit.times[1]){
        stop("tval is earlier than the first visit time")
    }
    
    #check if truncation.index is valid
    assert_that(is.count(truncation.index),truncation.index<=K)
    if(tval<=visit.times[truncation.index]){
        stop("tval is earlier than the left-truncation time")
    }
    
    K<-find.last.TRUE.index(visit.times<tval)
    if(length(models)!=K-truncation.index+1){
        stop("Number of models do not match number of time windows considered")
    }
    
    out<-list(covariate.data=covariate.data,formula=formula,visit.times=visit.times,tval=tval,truncation.index=truncation.index,models=models)
    class(out)<-"mult_stage_survfit"
    out
}


#' @export
print.mult_stage_survfit<-function(x,...){
    cat("\n")
    cat("mult_stage_survfit object\n")
    cat("Predictive model of survival probability with multiple stages\n\n")
    cat("visit times (visit.times): ",paste(x$visit.times,collapse=", "),"\n")
    cat("time for survival probability (tval): ",x$tval,"\n")
    cat("left-truncation index of visit times (truncation.index): ",x$truncation.index,"\n")
    cat("Formula: ",as.character(x$formula),"\n")
    cat(sprintf("survival probability models in time windows: list of %d model(s)",length(x$models)),"\n")
    if(inherits(x$models[[1]],"intercept_IF_model")){
        out<-predict(x)
        cat("Estimated survival probability for this intercept-only model: ",out$est,"\n")
        cat("Standard error of survival probability estimate: ",out$SE,"\n")
        cat(round(out$conf.level*100,2),"% confidence interval:\n")
        print(out$CI)
    }else if(inherits(x$models[[1]],"intercept_model")){
        cat("Estimated survival probability for this intercept-only model: ",predict(x),"\n")
    }
}


#' @title Predict and fitted method for mult_stage_survfit
#' @param object a "`mult_stage_survfit`" object
#' @param newdata new data frame of covariates used for prediction. Default is the `covariate.data` stored in `object`
#' @param pred.prob.trunc values at which to truncate predicted survival probabilities from each model in `object`. Must be a numeric vector with 2 numbers in \eqn{[0,1]}. Default is `c(0,1)`.
#' @param conf.level confidence level for confidence interval when using the SDR estimator to estimate a scalar estimand (rather than a function), i.e. a survival probability unconditional on covariates and possibly conditional on survival up to a visit time. Ignored otherwise. Default is 0.95.
#' @param ... further arguments passed to `predict` methods of the models in `object`
#' @return a vector of predicted (conditional) survival probabilities. If no covariates are used, that is, interested in the survival probability marginalized over all covariates so that all individuals have the same prediction, then only outputs that one predicted survial probability.
#' @export
predict.mult_stage_survfit<-function(object,newdata=object$covariate.data,pred.prob.trunc=c(0,1),conf.level=0.95,...){
    assert_that(is.numeric(pred.prob.trunc),length(pred.prob.trunc)==2)
    pred.prob.trunc<-sort(pred.prob.trunc)
    if(pred.prob.trunc[1]<0){
        stop("pred.prob.trunc must be >= 0")
    }
    if(pred.prob.trunc[2]>1){
        stop("pred.prob.trunc must be <= 1")
    }
    
    if(all(sapply(object$models,function(x) inherits(x,"intercept_IF_model")))){
        predictions<-sapply(object$models,function(model){
            pred<-model$value
            clip_interval(pred,pred.prob.trunc[1],pred.prob.trunc[2])
        })
        final.pred<-prod(predictions)
        
        #matrix to store IF in all time windows to facilitate computation of final IF. Each row is an observation, each col is a stage
        IF.mat<-matrix(0,nrow=nrow(object$covariate.data),ncol=length(object$models))
        rownames(IF.mat)<-names(object$models[[1]]$IF)
        
        for(k in 1:ncol(IF.mat)){
            IF<-object$models[[k]]$IF
            IF.mat[,k]<-IF[match(rownames(IF.mat),names(IF))]*nrow(IF.mat)/length(IF)
        }
        IF.mat<-ifelse(is.na(IF.mat),0,IF.mat)
        
        final.IF<-lapply(1:ncol(IF.mat),function(k){
            IF.mat[,k]*prod(predictions[-k])
        })%>%reduce(get("+"))
        
        SE<-sqrt(mean(final.IF^2))/sqrt(length(final.IF))
        qs<-c((1-conf.level)/2,(1+conf.level)/2)
        CI<-final.pred+qnorm(qs)*SE
        names(CI)<-paste(round(qs*100,1),"%")
        
        est_IF_SE_CI_conflevel(final.pred,final.IF,SE,CI,conf.level)
    }else{
        predictions<-lapply(object$models,function(model){
            if(inherits(model,"SuperLearner")){
                pred<-predict(model,newdata=newdata,...)$pred
            }else{
                pred<-predict(model,newdata=newdata,...)
            }
            clip_interval(pred,pred.prob.trunc[1],pred.prob.trunc[2])
        })
        reduce(predictions,get("*"))
    }
}

#' @rdname predict.mult_stage_survfit
#' @export
fitted.mult_stage_survfit<-function(object,pred.prob.trunc=c(0,1),...){
    predict.mult_stage_survfit(object=object,newdata=object$covariate.data,pred.prob.trunc=pred.prob.trunc,...)
}
