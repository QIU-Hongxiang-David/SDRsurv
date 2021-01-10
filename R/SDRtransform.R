#doubly robust transformation at one stage
#input fitted values should be truncated to minimal that can cover the time window to avoid unnecessary computation
#input tvals is sorted in increasing order; must all be greater than the smallest time in pred_event_censor_obj
#output: matrix of doubly transformed pseudo-outcome at the stage. each row is an individual; each column is a time in tvals
.DRtransform<-function(follow.up.time,pred_event_censor_obj,tvals,next.check.in.time=Inf,id.var,time.var,event.var,denom.survival.trunc){
    tvals.bar<-pmin(tvals,next.check.in.time) #tvals truncated at next check-in time
    
    #compute Ghat(s-) at s=event times
    Ghat.minus<-pred_event_censor_obj$event$surv
    for(i in 1:length(pred_event_censor_obj$event$time)){
        #j is the index of the first censoring time in censoring times that is >= current event time
        #will use j-1 to get Ghat(s-)
        #if no censoring time >= current event time, then set j to be the index after the last censoring time
        #if j=1 (all censoring time >= current event time), set Ghat(s-) to 1
        j<-find.first.TRUE.index(pred_event_censor_obj$censor$time>=pred_event_censor_obj$event$time[i],
                                 noTRUE=length(pred_event_censor_obj$censor$time)+1)
        if(j==1){
            Ghat.minus[,i]<-1
        }else{
            Ghat.minus[,i]<-pred_event_censor_obj$censor$surv[,j-1]
        }
    }
    
    #compute integrand in the DR transform at each event time
    integrand<-pred_event_censor_obj$event$surv
    for(i in 1:length(pred_event_censor_obj$event$time)){
        if(i==1){
            integrand[,i]<-(1-pred_event_censor_obj$event$surv[,i])/pred_event_censor_obj$event$surv[,i]/pmax(Ghat.minus[,i],denom.survival.trunc)
        }else{
            integrand[,i]<-(pred_event_censor_obj$event$surv[,i-1]-pred_event_censor_obj$event$surv[,i])/
                pred_event_censor_obj$event$surv[,i]/
                pred_event_censor_obj$event$surv[,i-1]/
                pmax(Ghat.minus[,i],denom.survival.trunc)
        }
    }
    #integral in the DR transform at each event time
    integral<-rowCumsums(integrand)
    
    output<-matrix(nrow=nrow(pred_event_censor_obj$event$surv),ncol=length(tvals.bar))
    rownames(output)<-rownames(pred_event_censor_obj$event$surv)
    colnames(output)<-as.character(tvals)
    for(i in 1:nrow(output)){
        id.matching.data<-follow.up.time%>%filter(.data[[id.var]]==rownames(output)[i])
        X<-pull(id.matching.data,.data[[time.var]])
        Delta<-pull(id.matching.data,.data[[event.var]])
        for(j in 1:ncol(output)){
            if(j>1 && tvals.bar[j]==tvals.bar[j-1]){
                output[i,j]<-output[i,j-1]
            }else{
                #find Shat at t
                #k.t is the index of the last event time in event times that is <= current tval (t)
                #will use k.t to get Shat at t
                #if no event time < t, then set Shat.t to be 1
                k.t<-find.last.TRUE.index(pred_event_censor_obj$event$time<=tvals.bar[j],noTRUE=0)
                if(k.t==0){
                    Shat.t<-1
                }else{
                    Shat.t<-pred_event_censor_obj$event$surv[i,k.t]
                }

                #find Shat at X
                #same logic as above
                k.SX<-find.last.TRUE.index(pred_event_censor_obj$event$time<=X,noTRUE=0)
                if(k.SX==0){
                    Shat.X<-1
                }else{
                    Shat.X<-pred_event_censor_obj$event$surv[i,k.SX]
                }
                
                #compute first IPW term in the bracket
                if(Delta==1 && X<=tvals.bar[j]){
                    #find Ghat(X-)
                    #same logic as above
                    k.GX<-find.first.TRUE.index(pred_event_censor_obj$censor$time>=X,
                                                noTRUE=length(pred_event_censor_obj$censor$time)+1)
                    if(k.GX==1){
                        Ghat.Xminus<-1
                    }else{
                        Ghat.Xminus<-pred_event_censor_obj$censor$surv[i,k.GX-1]
                    }
                    
                    IPW.term<-1/Shat.X/pmax(Ghat.Xminus,denom.survival.trunc)
                }else{
                    IPW.term<-0
                }
                
                #find integral at min(X,t)
                k.int<-min(k.t,k.SX)
                if(k.int==0){
                    integral.Xt<-0
                }else{
                    integral.Xt<-integral[i,k.int]
                }
                
                output[i,j]<-1-Shat.t
                if(Shat.t!=0){
                    output[i,j]<-output[i,j]+Shat.t*(IPW.term-integral.Xt)
                }
            }
        }
    }
    
    output
}


#' @title Sequentially doubly robust transformation of fitted survival and censoring probabilities
#' @name SDRtransform
#' @description Apply doubly robust transformation on fitted survival and censoring probabilities in each time window for estimating P(T <= t | T > truncation time, covariates available at truncation time).
#' @param follow.up.time see \code{\link{SDRsurv}}
#' @param pred_event_censor.list list of `pred_event_censor` objects. Each `pred_event_censor` object in the list corresponds to a time window in `check.in.times` after `truncation.index` in increasing order.
#' @param check.in.times see \code{\link{SDRsurv}}
#' @param tvals see \code{\link{SDRsurv}}. Must be sorted in ascending order.
#' @param truncation.index see \code{\link{SDRsurv}}
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param denom.survival.trunc see \code{\link{SDRsurv}}
#' @return  a list of list of data frames containing doubly robust transformed pseudo-outcomes that can be used by \code{\link{estQ.SuperLearner}}.
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{SDRsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
SDRtransform<-function(
    follow.up.time,
    pred_event_censor.list,
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
    
    index.shift<-truncation.index-1 #shift for the index of pred_event_censor.list
    
    #transformed is a list of matrices. Each matrix corresponds to a time window; each row is an individual followed up in the time window; each column is a time in tvals after start of time window
    transformed<-lapply(truncation.index:K,function(k){
        if(k<length(check.in.times)){
            pred_event_censor_obj<-truncate_pred_event_censor(pred_event_censor.list[[k-index.shift]],check.in.times[k+1])
        }else{
            pred_event_censor_obj<-pred_event_censor.list[[k-index.shift]]
        }
        
        if(k==K){
            next.check.in.time<-tail(tvals,1) #only need predictions up to the last tvals
        }else{
            next.check.in.time<-check.in.times[k+1]
        }
        valid.tvals<-tvals[tvals>check.in.times[k]]
        .DRtransform(follow.up.time,pred_event_censor_obj,valid.tvals,
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
