#' @title Estimate (conditional) survival probabilities with all three methods based on sequentially doubly robust, G-computation and inverse probability weighting (IPW) transformations
#' @name SDR_G_IPW_surv
#' @description
#' Estimate P(T > t | T > truncation time, covariates available at truncation time) for given t, where T is the time to event, using all three methods based on sequentially doubly robust, G-computation and inverse probability weighting (IPW) transformations. Use a user-specified flexible method to fit survival curves of time to event/censoring at each stage and then use \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} to regress pseudo-outcome on covariates in order to estimate P(T > t | T > truncation time, covariates available at truncation time). This function is faster than running \code{\link{SDRsurv}}, \code{\link{Gsurv}} and \code{\link{IPWsurv}} separately as the survival curves are fitted only once.
#'
#' @param covariates a list of data frames of covarates in the order of check in times. Each data frame contains the covariates collected at a check-in time. Data frames may have different numbers of variables (may collect different variables at different check-in times) and different numbers of individuals (some individuals may have an event or is censored before a later check-in time). All data frames must have a common character variable (see `id.var`) that identifies each individual but no other variables with common names. No missing data is allowed.
#' @param follow.up.time data frame of follow up times, i.e., times to event/censoring. Contains the variable that identifies each individual, the follow up times and an indicator of event/(right-)censoring. Follow up times must be numeric. Indicator of event/censoring should be binary with 0=censored, 1=event.
#' @param check.in.times numeric/integer vector of check-in times in ascending order. The first check-in time is typically the baseline.
#' @param tvals times t for which P(T > t) given covariates are computed (T is the time to event). Default is all unique event times in `follow.up.time`. Will be sorted in ascending order.
#' @param truncation.index index of the check-in time to which left-truncation is applied. The truncation time is `check.in.times[truncation.index]`. Covariates available up to (inclusive) `check.in.times[truncation.index]` are of interest. Default is 1, corresponding to no truncation.
#' @param id.var (character) name of the variable that identifies each individual.
#' @param time.var (character) name of the variable containing follow up times in the data frame `follow.up.time`.
#' @param event.var (character) name of the variable containing indicator of event/censoring in the data frame `follow.up.time`.
#' @param event.formula a list of formulas to specify covariates being used when estimating the conditional survival probabilities of time to event at each check-in time. The length should be the number of check in times after `truncation.index` (inclusive). Default is `~ .` for all check-in times, which includes main effects of all covariates available at each check-in time.
#' @param censor.formula a list of formulas to specify covariates being used when estimating the conditional survival probabilities of time to censoring at each check-in time. The length should be the number of check in times after `truncation.index` (inclusive). Default is `~ .` for all check-in times, which includes main effects of all covariates available at each check-in time.
#' @param Q.formula formula to specify covariates being used for estimating P(T > t | T > `check.in.times[truncation.index]`, covariates available at `check.in.times[truncation.index]`). Set to include intercept only (`~ 0` or `~ -1`) for marginal survival probability. Default is `~ .`, which includes main effects of all available covariates up to (inclusive) the `check.in.times[truncation.index]`.
#' @param event.method one of `"rfsrc"`, `"ctree"`, `"rpart"`, `"cforest"`, `"coxph"`, `"coxtime"`, `"deepsurv"`, `"dnnsurv"`, `"akritas"`. The machine learning method to fit survival survival curves of time to event in each time window. Default is `"rfsrc`. See the underlying wrappers \code{\link{fit_rfsrc}}, \code{\link{fit_ctree}}, \code{\link{fit_rpart}}, \code{\link{fit_cforest}}, \code{\link{fit_coxph}}, \code{\link{fit_coxtime}}, \code{\link{fit_deepsurv}}, \code{\link{fit_dnnsurv}}, \code{\link{fit_akritas}} for the available options.
#' @param censor.method one of `"rfsrc"`, `"ctree"`, `"rpart"`, `"cforest"`, `"coxph"`, `"coxtime"`, `"deepsurv"`, `"dnnsurv"`, `"akritas"`. The machine learning method to fit survival survival curves of time to censoring in each time window. Similar to `event.method`.
#' @param event.control a returned value from \code{\link{fit_surv_option}}
#' @param censor.control a returned value from \code{\link{fit_surv_option}}
#' @param Q.SuperLearner.control a list containing optional arguments passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}}. We encourage using a named list. Will be passed to \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} by running a command like `do.call(SuperLearner, Q.SuperLearner.control)`. Default is `list(SL.library="SL.lm")`, which uses linear regression. The user should not specify `Y` and `X`, and must specify `SL.library` if default is not used. If `family` is gaussian by default if unspecified, and must be gaussian if specified, with a possibly non-identity link. When `Q.formula` only includes an intercept, \code{\link[SuperLearner:SuperLearner]{SuperLearner::SuperLearner}} will not be called and the default setting can be used.
#' @param denom.survival.trunc the numeric truncation value for the denominator. All denominators below `denom.survival.trunc` will be set to `denom.survival.trunc` for numerical stability.
#' @return a list containing three lists of fitted `SuperLearner` models corresponding to each t in `tvals`. The three elements correspond to SDR, G-computation and IPW, respectively.
#' @section Formula arguments:
#' All formulas should have covariates on the right-hand side and no terms on the left-hand side, e.g., `~ V1 + V2 + V3`. At each check-in time, the corresponding formulas may (and usually should) contain covariates at previous check-in times, and must only include available covariates up to (inclusive) that check-in time. Interactions, polynomials and splines may be treated differently by different machine learning methods to estimate conditional survival curves.
#' @export
SDR_G_IPW_surv<-function(
    covariates,
    follow.up.time,
    check.in.times,
    tvals=NULL,
    truncation.index=1,
    id.var,
    time.var,
    event.var,
    event.formula=NULL,
    censor.formula=NULL,
    Q.formula=~.,
    event.method=c("rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","akritas"),
    censor.method=c("rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","akritas"),
    event.control=fit_surv_option(),
    censor.control=fit_surv_option(),
    Q.SuperLearner.control=list(family=gaussian(),SL.library="SL.lm"),
    denom.survival.trunc=1e-3
){
    assert_that(is.string(id.var))
    assert_that(is.string(time.var))
    assert_that(is.string(event.var))
    
    event.method<-match.arg(event.method)
    censor.method<-match.arg(censor.method)
    
    K<-length(check.in.times) #number of check-in times
    
    ############################################################################
    # check inputs are valid and set default values
    ############################################################################
    
    #check variables correctly exist
    if(!all(sapply(covariates,has_name,which=id.var))){
        stop(paste(id.var,"not present in 1+ covariates data"))
    }
    if(!all(sapply(covariates,function(d) is.character(pull(d,.data[[id.var]]))))){
        stop(paste(id.var,"is not character in 1+ covariates data"))
    }
    if(!has_name(follow.up.time,id.var)){
        stop(paste(id.var,"not present in follow.up.time"))
    }
    if(!is.character(pull(follow.up.time,.data[[id.var]]))){
        stop(paste(id.var,"is not character in follow.up.time"))
    }
    if(!has_name(follow.up.time,time.var)){
        stop(paste(time.var,"not present in follow.up.time"))
    }
    if(!has_name(follow.up.time,event.var)){
        stop(paste(event.var,"not present in follow.up.time"))
    }
    
    #check missing data
    if(!all(sapply(covariates,noNA))){
        stop("Missing data in 1+ covariates data")
    }
    if(!noNA(follow.up.time)){
        stop("Missing data in follow.up.time")
    }
    
    #check id.var is unique
    if(any(sapply(covariates,function(x) any(duplicated(pull(x,.data[[id.var]])))))){
        stop(paste("Duplicated",id.var,"in 1+ covariates data"))
    }
    if(any(duplicated(pull(follow.up.time,.data[[id.var]])))){
        stop(paste("Duplicated",id.var,"in follow.up.time"))
    }
    
    #check if time.var is numeric
    if(!is.numeric(pull(follow.up.time,.data[[time.var]]))){
        stop(paste(time.var),"is not numeric")
    }
    
    #check event.var is binary
    if(!(all(pull(follow.up.time,.data[[event.var]]) %in% c(0,1)))){
        stop(paste(event.var,"is not binary"))
    }
    
    all.event.times<-follow.up.time%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort

    #check if check.in.times are ascending with unique values
    if(is.unsorted(check.in.times,strictly=TRUE)){
        stop("check.in.times is not sorted in ascending order with unique values")
    }
    
    #check if all follow up times are >= the first check-in time
    if(!all(pull(follow.up.time,.data[[time.var]])>=check.in.times[1])){
        stop("At least one time in follow.up.time is earlier than the first check-in time")
    }
    
    #set default tvals and check they are numbers that are greater than the first check-in time
    if(is.null(tvals)){
        if(any(c(event.method,censor.method) %in% c("coxtime","deepsurv","dnnsurv","akritas"))){
            warning("When tvals are all event times, using coxtime, deepsurv, dnnsurv or akritas may lead to imprecision caused by conversion between numeric and character.")
        }
        tvals<-all.event.times
    }
    assert_that(is.numeric(tvals),noNA(tvals))
    tvals<-sort(tvals)
    if(!all(tvals>=check.in.times[1])){
        stop("At least one value in tvals is earlier than the first check-in time")
    }
    
    #check max(tvals) is reasonable
    if(tail(tvals,1)>max(all.event.times)){
        warning("At least one value in tvals is greater than the max time to event. Estimates on the tail may be non-informative.")
    }
    
    #check if truncation.index is valid
    assert_that(is.count(truncation.index),truncation.index<=K)
    if(!all(tvals>check.in.times[truncation.index])){
        stop("At least one value in tvals is earlier than the left-truncation time")
    }
    
    #check monotone missing of individuals
    if(K>1){
        lapply(2:K,function(i){
            if(!all(pull(covariates[[i]],.data[[id.var]]) %in% pull(covariates[[i-1]],.data[[id.var]]))){
                stop(paste0("1+ individual in covariates[[",i,"]] does not appear in covariates[[",i-1,"]]"))
            }
        })
    }
    
    #check individuals' follow up times are consistent with available covariates
    lapply(1:K,function(i){
        if(!setequal(pull(covariates[[i]],.data[[id.var]]),
                     follow.up.time%>%filter(.data[[time.var]]>check.in.times[i])%>%pull(.data[[id.var]]))){
            stop(paste0("Individuals in covariates[[",i,"]] differ from those being followed up after check.in.time[i]"))
        }
    })
    
    #check duplicate variable names in covaraites and follow.up.time
    lapply(1:K,function(i){
        if(length(intersect(setdiff(names(covariates[[i]]),id.var),
                            setdiff(names(follow.up.time),id.var)))>0){
            stop(paste0("Duplicated variables in covariates[[",i,"]] and follow.up.times"))
        }
    })
    
    #set default formulas for survival regressions and check if variables are all available at each check-in time
    #also check duplicated variable names in covariates
    if(is.null(event.formula)){
        event.formula<-lapply(check.in.times,function(x) ~.)
    }
    if(is.null(censor.formula)){
        censor.formula<-lapply(check.in.times,function(x) ~.)
    }
    history.covars<-NULL
    for(k in 1:K){
        if(any(names(covariates[[k]]) %in% history.covars)){
            stop("Duplicated variable names in covariates")
        }else{
            history.covars<-c(history.covars,setdiff(names(covariates[[k]]),id.var))
        }
        
        if(as.character(event.formula[[k]])[1]!="~"){
            stop(paste0("event.formula[[",k,"]] has variables on the left-hand side"))
        }
        if(as.character(censor.formula[[k]])[1]!="~"){
            stop(paste0("censor.formula[[",k,"]] has variables on the left-hand side"))
        }
        
        event.covars<-setdiff(all.vars(event.formula[[k]]),".")
        censor.covars<-setdiff(all.vars(censor.formula[[k]]),".")
        if(!all(event.covars %in% history.covars)){
            stop(paste0("event.formula[[",k,"]] contains varibales not available at check.in.times[",k,"]"))
        }
        if(!all(censor.covars %in% history.covars)){
            stop(paste0("censor.formula[[",k,"]] contains varibales not available at check.in.times[",k,"]"))
        }
    }
    
    #check if variables in Q.formula are all available at truncation time
    Q.covars<-setdiff(all.vars(Q.formula),".")
    history.covars<-setdiff(do.call(c,lapply(covariates[1:truncation.index],names)),id.var)
    if(!all(Q.covars %in% history.covars)){
        stop("Q.formula contains covariates not available up to truncation time")
    }
    
    #check if event.control and censor.control are fit_surv_option objects
    if(!inherits(event.control,"fit_surv_option")){
        stop("event.control is not a fit_surv_option object")
    }
    if(!inherits(censor.control,"fit_surv_option")){
        stop("censor.control is not a fit_surv_option object")
    }
    
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
    
    ############################################################################
    # run survival regressions
    ############################################################################
    #K is the last check.in.time that needs to be considered
    K<-find.last.TRUE.index(check.in.times<tail(tvals,1))
    
    index.shift<-truncation.index-1 #shift for the index of pred_event_censor.list
    
    pred_event_censor.list<-lapply(truncation.index:K,function(k){
        history<-reduce(covariates[1:k],.f=function(d1,d2){
            right_join(d1,d2,by=id.var)
        })%>%arrange(.data[[id.var]])
        
        if(k<length(check.in.times)){
            event.follow.up.time<-admin.censor(follow.up.time,time.var,event.var,check.in.times[k+1])
            censor.follow.up.time<-admin.censor(
                follow.up.time%>%
                    mutate("{event.var}":=1-.data[[event.var]],
                           "{time.var}":=left.shift.censoring(.data[[time.var]],.data[[event.var]])),
                time.var,event.var,check.in.times[k+1])
        }else{
            event.follow.up.time<-follow.up.time
            censor.follow.up.time<-follow.up.time%>%
                mutate("{event.var}":=1-.data[[event.var]],
                       "{time.var}":=left.shift.censoring(.data[[time.var]],.data[[event.var]]))
        }
        
        #fit event
        event.surv.data<-left_join(history,event.follow.up.time,by=id.var)
        form<-as.formula(paste("Surv(",time.var,",",event.var,")",
                               paste(as.character(event.formula[[k-index.shift]]),collapse=""),
                               collapse=""))
        fit_surv_arg<-c(
            list(method=event.method,formula=form,data=event.surv.data,id.var=id.var,time.var=time.var,event.var=event.var),
            event.control
        )
        pred_event_obj<-do.call(fit_surv,fit_surv_arg)
        
        
        #fir censoring
        censor.surv.data<-left_join(history,censor.follow.up.time,by=id.var)
        form<-as.formula(paste("Surv(",time.var,",",event.var,")",
                               paste(as.character(censor.formula[[k-index.shift]]),collapse=""),
                               collapse=""))
        fit_surv_arg<-c(
            list(method=censor.method,formula=form,data=censor.surv.data,id.var=id.var,time.var=time.var,event.var=event.var),
            censor.control
        )
        pred_censor_obj<-do.call(fit_surv,fit_surv_arg)
        
        
        pred_event_censor(pred_event_obj,pred_censor_obj)
    })
    
    ############################################################################
    # transformations and regression
    ############################################################################
    SDR.output<-SDRreg.SuperLearner(covariates,follow.up.time,pred_event_censor.list,check.in.times,tvals,truncation.index,id.var,time.var,event.var,Q.formula,Q.SuperLearner.control,denom.survival.trunc)
    G.output<-Greg.SuperLearner(covariates,follow.up.time,lapply(pred_event_censor.list,function(x) x$event),check.in.times,tvals,truncation.index,id.var,time.var,event.var,Q.formula,Q.SuperLearner.control)
    IPW.output<-IPWreg.SuperLearner(covariates,follow.up.time,lapply(pred_event_censor.list,function(x) x$censor),check.in.times,tvals,truncation.index,id.var,time.var,event.var,Q.formula,Q.SuperLearner.control,denom.survival.trunc)
    list(SDR=SDR.output,G=G.output,IPW=IPW.output)
}
