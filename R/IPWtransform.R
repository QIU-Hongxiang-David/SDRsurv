#IPW-computation transformation at one stage
#input fitted values should be truncated to minimal that can cover the time window to avoid unnecessary computation
#input tvals is sorted in increasing order; must all be greater than the smallest time in pred_censor_obj
#output: matrix of doubly transformed pseudo-outcome at the stage. each row is an individual; each column is a time in tvals
.IPWtransform<-function(follow.up.time,pred_censor_obj,tvals,next.check.in.time=Inf,id.var,time.var,event.var,denom.survival.trunc){
    tvals.bar<-pmin(tvals,next.check.in.time) #tvals truncated at next check-in time
    
    output<-matrix(nrow=nrow(pred_censor_obj$surv),ncol=length(tvals.bar))
    rownames(output)<-rownames(pred_censor_obj$surv)
    colnames(output)<-as.character(tvals)
    
    for(i in 1:nrow(output)){
        id.matching.data<-follow.up.time%>%filter(.data[[id.var]]==rownames(output)[i])
        X<-pull(id.matching.data,.data[[time.var]])
        Delta<-pull(id.matching.data,.data[[event.var]])
        for(j in 1:ncol(output)){
            if(j>1 && tvals.bar[j]==tvals.bar[j-1]){
                output[i,j]<-output[i,j-1]
            }else{
                if(Delta==1 && X<=tvals.bar[j]){
                    #find Ghat(X-)
                    #k.GX is the index of the first censoring time in censoring times that is >= X
                    #will use k.GX-1 to get Ghat(X-)
                    #if no censoring time >= X, then set k.GX to be the index after the last censoring time
                    #if k.GX=1 (all censoring time >= X), set Ghat(X-) to 1
                    k.GX<-find.first.TRUE.index(pred_censor_obj$time>=X,
                                                noTRUE=length(pred_censor_obj$time)+1)
                    if(k.GX==1){
                        Ghat.Xminus<-1
                    }else{
                        Ghat.Xminus<-pred_censor_obj$surv[i,k.GX-1]
                    }
                    
                    output[i,j]<-1/pmax(Ghat.Xminus,denom.survival.trunc)
                }else{
                    output[i,j]<-0
                }
            }
        }
    }
    
    output
}


#' @title Inverse probability weighting (IPW) transformation of fitted survival and censoring probabilities
#' @name IPWtransform
#' @description Apply inverse probability weighting (IPW) transformation on fitted survival and censoring probabilities in each time window for estimating P(T <= t | T > truncation time, covariates available at truncation time).
#' @param follow.up.time see \code{\link{SDRsurv}}
#' @param pred_censor.list list of `pred_surv` objects (see \code{\link{pred_surv}}) containing survival probabilities for time to censoring. Each `pred_surv` object in the list corresponds to a time window in `check.in.times` after `truncation.index` in increasing order.
#' @param check.in.times see \code{\link{SDRsurv}}
#' @param tvals see \code{\link{SDRsurv}}. Must be sorted in ascending order.
#' @param truncation.index see \code{\link{SDRsurv}}
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param denom.survival.trunc see \code{\link{SDRsurv}}
#' @return  a list of list of data frames containing doubly robust transformed pseudo-outcomes that can be used by \code{\link{estQ.SuperLearner}}.
#' @section Warning:
#' This function is designed to be called by \code{\link{IPWsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
IPWtransform<-function(
    follow.up.time,
    pred_censor.list,
    check.in.times,
    tvals,
    truncation.index,
    id.var,
    time.var,
    event.var,
    denom.survival.trunc=1e-3
){
    #K is the last check.in.time that needs to be considered
    K<-find.last.TRUE.index(check.in.times<tail(tvals,1))
    
    index.shift<-truncation.index-1 #shift for the index of pred_censor.list
    
    #transformed is a list of matrices. Each matrix corresponds to a time window; each row is an individual followed up in the time window; each column is a time in tvals after start of time window
    transformed<-lapply(truncation.index:K,function(k){
        if(k<length(check.in.times)){
            pred_censor_obj<-truncate_pred_surv(pred_censor.list[[k-index.shift]],check.in.times[k+1])
        }else{
            pred_censor_obj<-pred_censor.list[[k-index.shift]]
        }
        
        if(k==K){
            next.check.in.time<-tail(tvals,1) #only need predictions up to the last tvals
        }else{
            next.check.in.time<-check.in.times[k+1]
        }
        valid.tvals<-tvals[tvals>check.in.times[k]]
        .IPWtransform(follow.up.time,pred_censor_obj,valid.tvals,
                     next.check.in.time=next.check.in.time,
                     id.var,time.var,event.var,denom.survival.trunc)
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
