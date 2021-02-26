#' @title Inverse-probability weighting (IPW) transformation
#' @name IPWtransform
#' @description
#' Given a `pred_surv` object for time to censoring in a time window, calculates the IPW transformation in the time window. The transformation is used as the outcome when estimating the conditional survival probability at the next check-in time.
#' @param follow.up.time see \code{\link{SDRsurv}}
#' @param pred_censor_obj a `pred_surv` object for time to censoring in the time window of interest
#' @param tvals see \code{\link{SDRsurv}}. Must be sorted in ascending order and all greater than the smallest time in `pred_censor_obj`
#' @param next.check.in.time the next check-in time. Default is `Inf`, corresponding to the last time window
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param denom.survival.trunc see \code{\link{SDRsurv}}
#' @return a named matrix of transformations. Each row corresponds to an individual; each column corresponds to a value of `tvals`. Row names are elements in `follow.up.time$id.var`; column names are values of `tvals`.
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{IPWsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
IPWtransform<-function(follow.up.time,pred_censor_obj,tvals,next.check.in.time=Inf,id.var,time.var,event.var,denom.survival.trunc){
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
                    
                    output[i,j]<-1-1/pmax(Ghat.Xminus,denom.survival.trunc)
                }else{
                    output[i,j]<-1
                }
            }
        }
    }
    
    output
}


#' @title Regression based on IPW transformation of fitted survival and censoring probabilities using SuperLearner
#' @name IPWreg.SuperLearner
#' @description Apply IPW transformation on fitted survival and censoring probabilities in each time window and estimate P(T > t | T > truncation time, covariates available at truncation time) with \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}.
#' @param covariates see \code{\link{SDRsurv}}
#' @param follow.up.time see \code{\link{SDRsurv}}
#' @param pred_censor.list list of `pred_surv` objects for time to censoring. Each `pred_censor.list` object in the list corresponds to a time window in `check.in.times` after `truncation.index` in increasing order.
#' @param check.in.times see \code{\link{SDRsurv}}
#' @param tvals see \code{\link{SDRsurv}}. Must be sorted in ascending order.
#' @param truncation.index see \code{\link{SDRsurv}}
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param Q.formula formula to specify covariates being used for estimating P(T > t | T > `check.in.times[truncation.index]`, covariates available at `check.in.times[truncation.index]`). Set to include intercept only (`~ 0` or `~ -1`) for marginal survival probability, which is simply the mean of pseudo-outcomes. Default is `~ .`, which includes main effects of all available covariates up to (inclusive) the `truncation.time`.
#' @param Q.SuperLearner.control see \code{\link{SDRsurv}}
#' @param denom.survival.trunc see \code{\link{SDRsurv}}
#' @return  a list of \code{\link{mult_stage_survfit}} objects, each corresponding to a value in `tvals`
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{SDRsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @section Custom learners:
#' Custom learners may be specified by providing an element named `SL.library` in `Q.SuperLearner.control`.The user may refer to resources such as \url{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html} for a guide to create custom learners.
#' @export
IPWreg.SuperLearner<-function(
    covariates,
    follow.up.time,
    pred_censor.list,
    check.in.times,
    tvals,
    truncation.index,
    id.var,
    time.var,
    event.var,
    Q.formula=~.,
    Q.SuperLearner.control=list(family=gaussian(),SL.library="SL.lm"),
    denom.survival.trunc=1e-3
){
    assert_that(denom.survival.trunc>0,denom.survival.trunc<=1)
    
    #check if Q.SuperLearner.control is a list and whether it specifies Y or X
    assert_that(is.list(Q.SuperLearner.control))
    if(any(c("Y","X") %in% names(Q.SuperLearner.control))){
        stop("Q.SuperLearner.control should not not specify Y or X")
    }
    
    if(!("family" %in% names(Q.SuperLearner.control))){
        Q.SuperLearner.control$family<-gaussian()
    }
    if(is.character(Q.SuperLearner.control$family)){
        if(Q.SuperLearner.control$family!="gaussian"){
            stop("Q.SuperLearner.control$family must be gaussian")
        }
    }else if(is.function(Q.SuperLearner.control$family)){
        if(!all.equal(Q.SuperLearner.control$family,gaussian)){
            stop("Q.SuperLearner.control$family must be gaussian")
        }
    }else if(Q.SuperLearner.control$family$family!="gaussian"){
        stop("Q.SuperLearner.control$family must be gaussian")
    }
    
    if(!("SL.library" %in% names(Q.SuperLearner.control))){
        stop("Q.SuperLearner.control should specify SL.library")
    }
    
    
    #K is the last check.in.time that needs to be considered
    K<-find.last.TRUE.index(check.in.times<tail(tvals,1))
    
    index.shift<-truncation.index-1 #shift for the index of pred_censor.list
    
    history<-reduce(covariates[1:truncation.index],.f=function(d1,d2){
        right_join(d1,d2,by=id.var)
    })%>%arrange(.data[[id.var]])
    
    #stagewise.models is a list of list of models: the outer layer corresponds to time windows; the inner layer corresponds to tvals; a model is NULL if the corresponding t is before the time window
    stagewise.models<-lapply(truncation.index:K,function(k){
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
        
        models<-list()
        for(i in seq_along(tvals)){
            if(tvals[i]<=check.in.times[k]){ #no regression for this t
                models<-c(models,list(NULL))
            }else if(i>1 && tvals[i-1]>next.check.in.time){ #same time window of interest, same regression model as the previous t
                models<-c(models,models[i-1])
            }else{
                Y<-IPWtransform(follow.up.time,pred_censor_obj,tvals[i],
                                next.check.in.time=next.check.in.time,
                                id.var,time.var,event.var,denom.survival.trunc)
                X<-model.frame(Q.formula,data=history%>%filter(.data[[id.var]] %in% rownames(Y))%>%select(!.data[[id.var]]))
                Y<-Y[,1]
                if(ncol(X)==0){ #intercept-only model
                    models<-c(models,list(intercept_model(mean(Y))))
                }else{
                    SuperLearner.arg<-c(
                        list(Y=Y,X=X),
                        Q.SuperLearner.control
                    )
                    models<-c(models,list(do.call(SuperLearner,SuperLearner.arg)))
                }
            }
        }
        models
    })
    
    #rearrange stagewise.models to a list of mult_stage_survfit objects, each corresponding to a value in tvals
    lapply(seq_along(tvals),function(i){
        models<-lapply(stagewise.models,function(x){
            x[[i]]
        })
        models<-models[!sapply(models,is.null)]
        mult_stage_survfit(covariate.data=history,formula=Q.formula,check.in.times=check.in.times,tval=tvals[i],truncation.index=truncation.index,models=models)
    })
}
