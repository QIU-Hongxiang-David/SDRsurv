#' @title Constructor of S3 class `pred_surv`
#' @name pred_surv
#' @description Construct a `pred_surv` object containing predicted/fitted survival probabilities for multiple individuals in the same sample. Can only hold survival curves of the form of decreasing rightcontinuous step functions.
#' @param time collection of all unique event times where the predicted/fitted survival function may have a jump. Must be increasing.
#' @param surv matrix of predicted/fitted survival probabilities. Each row corresponds to an individual; each column corresponds to an event time in `time`. Row names should be `id.var` that identifies each individual (see \code{\link{SDRsurv}}).
#' @details If there is no event, set `time=Inf` and `surv` to be a matrix with only one column of 1.
#' @export
pred_surv<-function(time,surv){
    assert_that(is.numeric(time))
    if(is.unsorted(time,strictly=TRUE)){
        stop("time is not sorted in ascending order with unique values")
    }
    assert_that(length(time)==ncol(surv))
    if(is.null(rownames(surv))){
        warning("surv does not have row names. This may cause trouble when linking predicted/fitted values in later usage.")
    }
    
    output<-list(time=time,surv=surv)
    class(output)<-"pred_surv"
    output
}


#' @title Constructor of S3 class `pred_event_censor`
#' @name pred_event_censor
#' @description Construct a `pred_event_censor` object containing predicted/fitted survival and censoring probabilities for multiple individuals in the same sample.  Can only hold survival curves of the form of decreasing right-continuous step functions.
#' @param pred_event_obj a `pred_surv` object holding predicted/fitted survival probabilities for events.
#' @param pred_censor_obj a `pred_surv` object holding predicted/fitted survival probabilities for censoring. Set to `NULL` (default) if no censoring present.
#' @details The predicted/fitted probability matrices in `pred_event_obj` and `pred_censor_obj` must have same row names, i.e., same number of individuals.
#' @export
pred_event_censor<-function(pred_event_obj,pred_censor_obj=NULL){
    if(!("pred_surv" %in% class(pred_event_obj))){
        stop("pred_event_obj is not a pred_surv object")
    }
    
    if(is.null(pred_censor_obj)){
        mat<-matrix(1,nrow=length(pred_event_obj$time),ncol=1)
        rownames(mat)<-rownames(pred_event_obj$surv)
        pred_censor_obj<-pred_surv(Inf,mat)
    }else if(!("pred_surv" %in% class(pred_censor_obj))){
        stop("pred_cens_obj is not a pred_surv object")
    }
    
    if(!setequal(rownames(pred_event_obj$surv),rownames(pred_censor_obj$surv))){
        stop("Different individuals in pred_event_obj and pred_censor_obj")
    }
    
    output<-list(event=pred_event_obj,censor=pred_censor_obj)
    class(output)<-"pred_event_censor"
    output
}


#' @title Truncate `pred_surv` object to minimal that covers a time window
#' @name truncate_pred_surv
#' @description Truncate the time and survival probabilities in a `pred_surv` object to minimal so that predicted/fitted survival probability at `end` is kept. This function can be used before applying transformations (e.g., \code{\link{SDRtransform}}, \code{\link{Gtransform}}, \code{\link{IPWtransform}}) to avoid unnecessary computation. No truncation is applied if `end` is greater than or equal to the last time in the `pred_surv` object.
#' @param pred_surv_obj `pred_surv` object.
#' @param end ending point (inclusive) of at which truncation applies. Default is `Inf`.
#' @return A `pred_surv` object with truncated times.
#' @export
truncate_pred_surv<-function(pred_surv_obj,end=Inf){
    #i is the index of the first time that is >= end
    #if none >= end, return the original object
    i<-find.first.TRUE.index(pred_surv_obj$time>=end,noTRUE=length(pred_surv_obj$time))
    
    if(i==length(pred_surv_obj$time)){
        pred_surv_obj
    }else{
        pred_surv(pred_surv_obj$time[1:i],pred_surv_obj$surv[,1:i,drop=FALSE])
    }
}

#' @title Truncate `pred_event_censor` object to minimal that covers a time window
#' @name truncate_pred_event_censor
#' @description Truncate the time and survival probabilities in a `pred_event_censor` object to minimal so that predicted/fitted survival and censoring probabilities at `end` are kept. This function can be used before applying transformations (e.g., \code{\link{SDRtransform}}, \code{\link{Gtransform}}, \code{\link{IPWtransform}}) to avoid unnecessary computation. No truncation is applied if `end` is greater than or equal to the last time in the `pred_event_censor` object.
#' @param pred_event_censor_obj `pred_event_censor` object.
#' @param end ending point (inclusive) of at which truncation applies. Default is `Inf`.
#' @return A `pred_event_censor` object with truncated times.
#' @export
truncate_pred_event_censor<-function(pred_event_censor_obj,end=Inf){
    pred_event_censor(truncate_pred_surv(pred_event_censor_obj$event,end),
                   truncate_pred_surv(pred_event_censor_obj$censor,end))
}
