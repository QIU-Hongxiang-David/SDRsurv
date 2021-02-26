#' @title S3 class for predictive models with an intercept only
#' @name intercept_model
#' @param value: the predicted mean value
#' @return an "`intercept_model`" object, essentially a list with element `value`.
#' @export
intercept_model<-function(value){
    assert_that(is.number(value))
    out<-list(value=value)
    class(out)<-"intercept_model"
    out
}

predict.intercept_model<-function(object,...){
    object$value
}

fitted.intercept_model<-predict.intercept_model


#' @title S3 class for the predictive model for (conditional) survival probability with multiple stages
#' @name mult_stage_survfit
#' @description
#' This function is the creator of `mult_stage_survfit` class. A `mult_stage_survfit` object is a collection of fitted models for conditional survival probability within each time window. The associated `predict` and `fitted` methods autocmatically conduct the prediction, which is multiplying all model predictions, potentially after clipping the predictions into an interval in \eqn{[0,1]}.
#' @param covariate.data data frame containing the covariates used to train the models
#' @param formula formula used in the models
#' @param check.in.times numeric/integer vector of check-in times in ascending order. The first check-in time is typically the baseline.
#' @param tvals time t for which P(T > t) (potentially given covariates) is estimated (T is the time to event)
#' @param truncation.index index of the check-in time to which left-truncation is applied. The truncation time is `check.in.times[truncation.index]`. Covariates available up to (inclusive) `check.in.times[truncation.index]` are of interest.
#' @param models list of models. Each model can predict the survival probability beyond the next check-in time given survival at the current check-in time. Models should be in an ascending order according to the corresponding check-in times.
#' @return a "`mult_stage_survfit`" object, essentially a list with input parameters of this function.
#' @export
mult_stage_survfit<-function(covariate.data,formula,check.in.times,tval,truncation.index,models){
    K<-length(check.in.times)
    
    assert_that(is.list(models))
    assert_that(is_formula(formula))
    
    #check if check.in.times are ascending with unique values
    if(is.unsorted(check.in.times,strictly=TRUE)){
        stop("check.in.times is not sorted in ascending order with unique values")
    }
    
    assert_that(is.number(tval),noNA(tval))
    if(tval<check.in.times[1]){
        stop("tval is earlier than the first check-in time")
    }
    
    #check if truncation.index is valid
    assert_that(is.count(truncation.index),truncation.index<=K)
    if(tval<=check.in.times[truncation.index]){
        stop("tval is earlier than the left-truncation time")
    }
    
    K<-find.last.TRUE.index(check.in.times<tval)
    if(length(models)!=K-truncation.index+1){
        stop("Number of models do not match number of time windows considered")
    }
    
    out<-list(covariate.data=covariate.data,formula=formula,check.in.times=check.in.times,tval=tval,truncation.index=truncation.index,models=models)
    class(out)<-"mult_stage_survfit"
    out
}

#' @export
print.mult_stage_survfit<-function(x,...){
    cat("\n")
    cat("mult_stage_survfit object\n")
    cat("Predictive model of survival probability with multiple stages\n\n")
    cat("check-in times (check.in.times): ",paste(x$check.in.times,collapse=", "),"\n")
    cat("time for survival probability (tval): ",x$tval,"\n")
    cat("left-truncation index of check-in times (truncation.index): ",x$truncation.index,"\n")
    cat("Formula: ",as.character(x$formula),"\n")
    cat(sprintf("survival probability models in time windows: list of %d model(s)",length(x$models)),"\n")
    if(inherits(x$models[[1]],"intercept_model")){
        cat("Estimated survival probability for this intercept-only model: ",predict(x),"\n")
    }
}


#' @title Predict and fitted method for mult_stage_survfit
#' @param object a "`mult_stage_survfit`" object
#' @param newdata new data frame of covariates used for prediction. Default is the `covariate.data` stored in `object`
#' @param pred.prob.trunc values at which to truncate predicted survival probabilities from each model in `object`. Must be a numeric vector with 2 numbers in \eqn{[0,1]}. Default is `c(0,1)`.
#' @param ... further arguments passed to `predict` methods of the models in `object`
#' @return a vector of predicted (conditional) survival probabilities. If no covariates are used, that is, interested in the survival probability marginalized over all covariates so that all individuals have the same prediction, then only outputs that one predicted survial probability.
#' @export
predict.mult_stage_survfit<-function(object,newdata=object$covariate.data,pred.prob.trunc=c(0,1),...){
    assert_that(is.numeric(pred.prob.trunc),length(pred.prob.trunc)==2)
    pred.prob.trunc<-sort(pred.prob.trunc)
    if(pred.prob.trunc[1]<0){
        stop("pred.prob.trunc must be >= 0")
    }
    if(pred.prob.trunc[2]>1){
        stop("pred.prob.trunc must be <= 1")
    }
    
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

#' @rdname predict.mult_stage_survfit
#' @export
fitted.mult_stage_survfit<-function(object,pred.prob.trunc=c(0,1),...){
    predict.mult_stage_survfit(object=object,newdata=object$covariate.data,pred.prob.trunc=pred.prob.trunc,...)
}
