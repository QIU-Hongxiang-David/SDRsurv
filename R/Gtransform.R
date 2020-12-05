#G-computation transformation at one stage
#input fitted values should be truncated to minimal that can cover the time window to avoid unnecessary computation
#input tvals is sorted in increasing order; must all be greater than the smallest time in pred_event_obj
#output: matrix of doubly transformed pseudo-outcome at the stage. each row is an individual; each column is a time in tvals
.Gtransform<-function(follow.up.time,pred_event_obj,tvals,next.check.in.time=Inf,id.var,time.var,event.var){
    tvals.bar<-pmin(tvals,next.check.in.time) #tvals truncated at next check-in time
    
    output<-matrix(nrow=nrow(pred_event_obj$surv),ncol=length(tvals.bar))
    rownames(output)<-rownames(pred_event_obj$surv)
    colnames(output)<-as.character(tvals)
    
    for(j in 1:ncol(output)){
        if(j>1 && tvals.bar[j]==tvals.bar[j-1]){
            output[,j]<-output[,j-1]
        }else{
            #find Shat at t
            #k.t is the index of the last event time in event times that is <= current tval (t)
            #will use k.t to get Shat at t
            #if no event time < t, then set Shat.t to be 1
            k.t<-find.last.TRUE.index(pred_event_obj$time<=tvals.bar[j],noTRUE=0)
            if(k.t==0){
                Shat<-1
            }else{
                Shat<-pred_event_obj$surv[,k.t]
            }
            output[,j]<-1-Shat
        }
    }
    output
}


#' @title G-computation transformation of fitted survival and censoring probabilities
#' @name Gtransform
#' @description Apply G-computation transformation on fitted survival and censoring probabilities in each time window for estimating P(T <= t | T > truncation time, covariates available at truncation time).
#' @param follow.up.time see \code{\link{SDRsurv}}
#' @param pred_event.list list of `pred_surv` objects (see \code{\link{pred_surv}}) containing survival probabilities for time to event. Each `pred_surv` object in the list corresponds to a time window in `check.in.times` after `truncation.index` in increasing order.
#' @param check.in.times see \code{\link{SDRsurv}}
#' @param tvals see \code{\link{SDRsurv}}. Must be sorted in ascending order.
#' @param truncation.index see \code{\link{SDRsurv}}
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @return  a list of list of data frames containing doubly robust transformed pseudo-outcomes that can be used by \code{\link{estQ.SuperLearner}}.
#' @section Warning:
#' This function is designed to be called by \code{\link{Gsurv}}. Therefore, inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
Gtransform<-function(
    follow.up.time,
    pred_event.list,
    check.in.times,
    tvals,
    truncation.index,
    id.var,
    time.var,
    event.var
){
    #K is the last check.in.time that needs to be considered
    K<-find.last.TRUE.index(check.in.times<tail(tvals,1))
    
    index.shift<-truncation.index-1 #shift for the index of pred_event.list
    
    #transformed is a list of matrices. Each matrix corresponds to a time window; each row is an individual followed up in the time window; each column is a time in tvals after start of time window
    transformed<-lapply(truncation.index:K,function(k){
        if(k<length(check.in.times)){
            pred_event_obj<-truncate_pred_surv(pred_event.list[[k-index.shift]],check.in.times[k+1])
        }else{
            pred_event_obj<-pred_event.list[[k-index.shift]]
        }
        
        if(k==K){
            next.check.in.time<-tail(tvals,1) #only need predictions up to the last tvals
        }else{
            next.check.in.time<-check.in.times[k+1]
        }
        valid.tvals<-tvals[tvals>check.in.times[k]]
        .Gtransform(follow.up.time,pred_event_obj,valid.tvals,
                     next.check.in.time=next.check.in.time,
                     id.var,time.var,event.var)
    })
    
    stagewise.pseudo.outcomes<-lapply(tvals,function(t){
        #J is the last check.in.time that needs to be considered
        J<-find.last.TRUE.index(check.in.times<t)
        lapply(truncation.index:J,function(j){
            tibble("{id.var}":=rownames(transformed[[j-index.shift]]),
                   pseudo.outcome=transformed[[j-index.shift]][,as.character(t)])
        })
    })
    names(stagewise.pseudo.outcomes)<-as.character(tvals)
    stagewise.pseudo.outcomes
}
