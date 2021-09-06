#' @title Options of machine learning methods' wrappers for fitting conditional survival curves
#' @name fit_surv_option
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to the wrapped machine learning function. Will be used in a command like `do.call(machine.learning, option)` where `machine.learning` is the machine learning function being called. `formula` and `data` should not be specified. For \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}, if `tune=TRUE`, then `mtry` and `nodesize` should not be specified either.
#' @param oob whether to use out-of-bag (OOB) fitted values from random forests (\code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} and \code{\link[party:cforest]{party::cforest}}) when sample splitting is not used (`nfold=1`). Ignored otherwise.
#' @param tune whether to tune `mtry` and `nodesize` for \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}. Ignored for other methods.
#' @param tune.option a list containing optional arguments passed to \code{\link[randomForestSRC:tune]{randomForestSRC::tune.rfsrc}} if \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} is used and `tune=TRUE`; ignored otherwise. `doBest` should not be specified.
#' @param lambda bandwidth parameter for uniform smoothing kernel in nearest neighbours estimation for method `"akritas"`. The default value of 0.5 is arbitrary and should be chosen by the user
#' @export
fit_surv_option<-function(nfold=1,option=list(),oob=TRUE,tune=TRUE,tune.option=list(),lambda=0.5){
    assert_that(is.count(nfold))
    assert_that(is.flag(oob))
    assert_that(is.flag(tune))
    assert_that(is.number(lambda),lambda>0)
    out<-list(nfold=nfold,option=option,oob=oob,tune=tune,tune.option=tune.option,lambda=lambda)
    class(out)<-"fit_surv_option"
    out
}


fit_surv<-function(method=c("survSuperLearner","rfsrc","ctree","rpart","cforest","coxph","coxtime","deepsurv","dnnsurv","akritas"),...){
    method<-match.arg(method)
    if(method=="survSuperLearner"){
        fit_survSuperLearner(...)
    }else if(method=="rfsrc"){
        fit_rfsrc(...)
    }else if(method=="ctree"){
        fit_ctree(...)
    }else if(method=="rpart"){
        fit_rpart(...)
    }else if(method=="cforest"){
        fit_cforest(...)
    }else if(method=="coxph"){
        fit_coxph(...)
    }else if(method=="coxtime"){
        fit_coxtime(...)
    }else if(method=="deepsurv"){
        fit_deepsurv(...)
    }else if(method=="dnnsurv"){
        fit_dnnsurv(...)
    }else if(method=="akritas"){
        fit_akritas(...)
    }
}


fit_no_event<-function(data,id.var){
    surv<-matrix(1,nrow=nrow(data),ncol=1)
    rownames(surv)<-pull(data,.data[[id.var]])
    pred_surv(Inf,surv)
}


#' @title Wrapper of `survSuperLearner::survSuperLearner`
#' @name fit_survSuperLearner
#' @param formula formula containing all covariates to be used
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survSuperLearner:survSuperLearner]{survSuperLearner::survSuperLearner}}. We encourage using a named list. Will be passed to \code{\link[survSuperLearner:survSuperLearner]{survSuperLearner::survSuperLearner}} by running a command like `do.call(survSuperLearner, option)`. The user should not specify `time`, `event`, `X`, `newX` or `new.times`. We encourage the user to specify `event.SL.library` and `cens.SL.library`.
#' @param ... ignored
#' @return a \code{\link{pred_event_censor}} class containing fitted survival curves for individuals in `data`
#' @export
fit_survSuperLearner<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(event.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc"),cens.SL.library=c("survSL.coxph","survSL.weibreg","survSL.gam","survSL.rfsrc")),...){
    .requireNamespace("survSuperLearner")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("time","event","X","newX","new.times") %in% names(option))){
        stop("option specifies time, event, X, newX or new.times")
    }
    
    if(nfold==1){
        time<-data%>%pull(time.var)
        event<-data%>%pull(event.var)
        newX<-X<-model.frame(formula,data=data%>%select(!c(.data[[id.var]],.data[[time.var]],.data[[event.var]])))
        new.times<-sort(unique(time))
        
        arg<-c(list(time=time,event=event,X=X,newX=newX,new.times=new.times),option)
        model<-do.call(survSuperLearner::survSuperLearner,arg)
        
        event.pred<-model$event.SL.predict
        row.names(event.pred)<-data%>%pull(id.var)
        
        censor.pred<-model$cens.SL.predict
        row.names(censor.pred)<-data%>%pull(id.var)
        pred_event_censor(pred_surv(time=new.times,surv=event.pred),
                          pred_surv(time=new.times,surv=censor.pred))
    }else{
        all.times<-data%>%pull(.data[[time.var]])%>%unique%>%sort
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        
        pred_event_censor.list<-lapply(folds,function(fold){
            d<-data%>%filter(!(.data[[id.var]] %in% .env$fold))
            test.d<-data%>%filter(.data[[id.var]] %in% .env$fold)
            
            time<-d%>%pull(time.var)
            event<-d%>%pull(event.var)
            X<-model.frame(formula,data=d%>%select(!c(.data[[id.var]],.data[[time.var]],.data[[event.var]])))
            newX<-model.frame(formula,data=test.d%>%select(!c(.data[[id.var]],.data[[time.var]],.data[[event.var]])))
            new.times<-all.times
            
            arg<-c(list(time=time,event=event,X=X,newX=newX,new.times=new.times),option)
            model<-do.call(survSuperLearner::survSuperLearner,arg)
            
            event.pred<-model$event.SL.predict
            row.names(event.pred)<-fold
            
            censor.pred<-model$cens.SL.predict
            row.names(censor.pred)<-fold
            
            pred_event_censor(pred_surv(all.times,event.pred),pred_surv(all.times,censor.pred))
        })
        
        event.pred<-lapply(pred_event_censor.list,function(x){
            x$event$surv
        })%>%do.call(what=rbind)
        event.pred<-event.pred[order(rownames(event.pred)),,drop=FALSE]
        
        censor.pred<-lapply(pred_event_censor.list,function(x){
            x$censor$surv
        })%>%do.call(what=rbind)
        censor.pred<-censor.pred[order(rownames(censor.pred)),,drop=FALSE]
        
        pred_event_censor(pred_surv(all.times,event.pred),pred_surv(all.times,censor.pred))
    }
}


#' @title Wrapper of `randomForestSRC::rfsrc`
#' @name fit_rfsrc
#' @param formula formula used by \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}. We encourage using a named list. Will be passed to \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} by running a command like `do.call(rfsrc, option)`. The user should not specify `formula` and `data`.
#' @param oob whether to use out-of-bag (OOB) fitted values from \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} when sample splitting is not used (`nfold=1`)
#' @param tune whether to tune `mtry` and `nodesize`.
#' @param tune.option a list containing optional arguments passed to \code{\link[randomForestSRC:tune]{randomForestSRC::tune.rfsrc}} if `tune=TRUE`; ignored otherwise. `doBest` should not be specified.
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_rfsrc<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),oob=TRUE,tune=TRUE,tune.option=list(),...){
    .requireNamespace("randomForestSRC")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data") %in% names(option))){
        stop("option specifies formula or data")
    }
    
    #check if option is a list and whether it specifies doBest
    if(tune){
        assert_that(is.list(tune.option))
        if("doBest" %in% names(tune.option)){
            stop("tune.option specifies doBest")
        }
        if(any(c("mtry","nodesize") %in% names(option))){
            stop("option specifies mtry or nodesize with tune=TRUE")
        }
    }
    
    #check if oob is logical
    assert_that(is.flag(oob))
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        if(tune){
            tune.arg<-c(
                list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
                tune.option
            )
            tune.output<-do.call(randomForestSRC::tune.rfsrc,tune.arg)
        }
        
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        if(tune){
            arg<-c(arg,list(mtry=tune.output$optimal["mtry"],nodesize=tune.output$optimal["nodesize"]))
        }
        model<-do.call(randomForestSRC::rfsrc,arg)
        if(oob){
            surv<-model$survival.oob
        }else{
            surv<-model$survival
        }
        rownames(surv)<-pull(data,.data[[id.var]])
        pred_surv(time=model$time.interest,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                if(tune){
                    tune.arg<-c(
                        list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                        tune.option
                    )
                    tune.output<-do.call(randomForestSRC::tune.rfsrc,tune.arg)
                }
                
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                if(tune){
                    arg<-c(arg,list(mtry=tune.output$optimal["mtry"],nodesize=tune.output$optimal["nodesize"]))
                }
                model<-do.call(randomForestSRC::rfsrc,arg)
                predict.model<-predict(model,data%>%filter(.data[[id.var]] %in% .env$fold))
                surv<-lapply(all.times,function(t){
                    i<-find.first.TRUE.index(predict.model$time.interest<=t,noTRUE=0)
                    if(i==0){
                        out<-matrix(1,nrow=length(fold),ncol=1)
                    }else{
                        out<-predict.model$survival[,i]
                    }
                    out
                })%>%do.call(what=cbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}




#' @title Wrapper of `party::ctree`
#' @name fit_ctree
#' @param formula formula used by \code{\link[party:ctree]{party::ctree}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[party:ctree]{party::ctree}}. We encourage using a named list. Will be passed to \code{\link[party:ctree]{party::ctree}} by running a command like `do.call(ctree, option)`. The user should not specify `formula` and `data`.
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_ctree<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),...){
    .requireNamespace("party")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data") %in% names(option))){
        stop("option specifies formula or data")
    }
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        model<-do.call(party::ctree,arg)
        surv<-lapply(party::treeresponse(model),function(surv_fit){
            sapply(all.times,function(t){
                i<-find.first.TRUE.index(surv_fit$time<=t,noTRUE=0)
                if(i==0){
                    out<-1
                }else{
                    out<-surv_fit$surv[i]
                }
                out
            })
        })%>%do.call(what=rbind)
        rownames(surv)<-pull(data,.data[[id.var]])
        pred_surv(time=all.times,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                model<-do.call(party::ctree,arg)
                surv<-lapply(party::treeresponse(model,newdata=data%>%filter(.data[[id.var]] %in% .env$fold)),function(surv_fit){
                    sapply(all.times,function(t){
                        i<-find.first.TRUE.index(surv_fit$time<=t,noTRUE=0)
                        if(i==0){
                            out<-1
                        }else{
                            out<-surv_fit$surv[i]
                        }
                        out
                    })
                })%>%do.call(what=rbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}



#' @title Wrapper of `rpart::rpart` and `partykit::as.party`
#' @name fit_rpart
#' @description First use \code{\link[rpart:rpart]{rpart::rpart}} to obtain a tree and then use \code{\link[partykit:party-coercion]{as.party}} to obtain Kaplan-Meier fits.
#' @param formula formula used by \code{\link[rpart:rpart]{rpart::rpart}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[rpart:rpart]{rpart::rpart}}. We encourage using a named list. Will be passed to \code{\link[rpart:rpart]{rpart::rpart}} by running a command like `do.call(rpart, option)`. The user should not specify `formula` and `data`.
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_rpart<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),...){
    .requireNamespace("rpart")
    .requireNamespace("party")
    .requireNamespace("partykit")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data") %in% names(option))){
        stop("option specifies formula or data")
    }
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        model<-partykit::as.party(do.call(rpart::rpart,arg))
        surv<-lapply(predict(model,type="prob"),function(surv_fit){
            sapply(all.times,function(t){
                i<-find.first.TRUE.index(surv_fit$time<=t,noTRUE=0)
                if(i==0){
                    out<-1
                }else{
                    out<-surv_fit$surv[i]
                }
                out
            })
        })%>%do.call(what=rbind)
        rownames(surv)<-pull(data,.data[[id.var]])
        pred_surv(time=all.times,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                model<-partykit::as.party(do.call(rpart::rpart,arg))
                surv<-lapply(predict(model,newdata=data%>%filter(.data[[id.var]] %in% .env$fold),type="prob"),function(surv_fit){
                    sapply(all.times,function(t){
                        i<-find.first.TRUE.index(surv_fit$time<=t,noTRUE=0)
                        if(i==0){
                            out<-1
                        }else{
                            out<-surv_fit$surv[i]
                        }
                        out
                    })
                })%>%do.call(what=rbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}




#' @title Wrapper of `party::cforest`
#' @name fit_cforest
#' @param formula formula used by \code{\link[party:cforest]{party::cforest}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[party:cforest]{party::cforest}}. We encourage using a named list. Will be passed to \code{\link[party:cforest]{party::cforest}} by running a command like `do.call(cforest, option)`. The user should not specify `formula` and `data`.
#' @param oob whether to use out-of-bag (OOB) fitted values from \code{\link[party:cforest]{party::cforest}} when sample splitting is not used (`nfold=1`)
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_cforest<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),oob=TRUE,...){
    .requireNamespace("party")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data") %in% names(option))){
        stop("option specifies formula or data")
    }
    
    #check if oob is logical
    assert_that(is.flag(oob))
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        model<-do.call(party::cforest,arg)
        surv<-lapply(party::treeresponse(model,OOB=oob),function(surv_fit){
            sapply(all.times,function(t){
                i<-find.first.TRUE.index(surv_fit$time<=t,noTRUE=0)
                if(i==0){
                    out<-1
                }else{
                    out<-surv_fit$surv[i]
                }
                out
            })
        })%>%do.call(what=rbind)
        rownames(surv)<-pull(data,.data[[id.var]])
        pred_surv(time=all.times,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                model<-do.call(party::cforest,arg)
                surv<-lapply(party::treeresponse(model,newdata=data%>%filter(.data[[id.var]] %in% .env$fold)),function(surv_fit){
                    sapply(all.times,function(t){
                        i<-find.first.TRUE.index(surv_fit$time<=t,noTRUE=0)
                        if(i==0){
                            out<-1
                        }else{
                            out<-surv_fit$surv[i]
                        }
                        out
                    })
                })%>%do.call(what=rbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}




#' @title Wrapper of `survival::coxph`
#' @name fit_coxph
#' @param formula formula used by \code{\link[survival:coxph]{survival::coxph}}. Currently \code{\link[survival]{strata}} is not supported.
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survival:coxph]{survival::coxph}}. We encourage using a named list. Will be passed to \code{\link[survival:coxph]{survival::coxph}} by running a command like `do.call(coxph, option)`. The user should not specify `formula` and `data`.
#' @param ... ignored
#' @param option a list containing optional arguments passed to \code{\link[survival:coxph]{survival::coxph}}. We encourage using a named list. Will be passed to \code{\link[survival:coxph]{survival::coxph}} by running a command like `do.call(coxph, option)`. The user should not specify `formula` and `data`.
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_coxph<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),...){
    .require("survival")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data") %in% names(option))){
        stop("option specifies formula or data")
    }
    
    if(any(grepl("strata\\(.*\\)",as.character(formula)))){
        warning("strata() function seems to be used in the formula. This may lead to an error.")
    }
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        cox.model<-do.call(survival::coxph,arg)
        surv<-lapply(all.times,function(t){
            predict(cox.model,newdata=data%>%mutate("{time.var}":=t),type="survival")
        })%>%do.call(what=cbind)
        # model<-survival::survfit(cox.model,newdata=data)
        # surv<-t(as.matrix.rowvec(model$surv))
        rownames(surv)<-pull(data,.data[[id.var]])
        pred_surv(time=all.times,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                cox.model<-do.call(survival::coxph,arg)
                surv<-lapply(all.times,function(t){
                    predict(cox.model,newdata=data%>%filter(.data[[id.var]] %in% .env$fold)%>%mutate("{time.var}":=t),type="survival")
                })%>%do.call(what=cbind)
                # model<-survival::survfit(cox.model,newdata=data%>%filter(.data[[id.var]] %in% .env$fold))
                # model.surv<-t(as.matrix.rowvec(model$surv))
                # surv<-lapply(all.times,function(t){
                #     i<-find.first.TRUE.index(model$time.interest<=t,noTRUE=0)
                #     if(i==0){
                #         out<-matrix(1,nrow=length(fold),ncol=1)
                #     }else{
                #         out<-as.matrix(model$surv)[,i]
                #     }
                #     out
                # })%>%do.call(what=cbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}




#' @title Wrapper of `survivalmodels::coxtime`
#' @name fit_coxtime
#' @param formula formula used by \code{\link[survivalmodels:coxtime]{survivalmodels::coxtime}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survivalmodels:coxtime]{survivalmodels::coxtime}}. We encourage using a named list. Will be passed to \code{\link[survivalmodels:coxtime]{survivalmodels::coxtime}} by running a command like `do.call(coxtime, option)`. The user should not specify `formula`, `data` and `reverse`; `time_variable`, `status_variable`, `x`, `y` will be ignored.
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_coxtime<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),...){
    .requireNamespace("survivalmodels")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data","reverse") %in% names(option))){
        stop("option specifies formula, data or reverse")
    }
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        model<-do.call(survivalmodels::coxtime,arg)
        surv<-predict(model,type="survival",distr6=FALSE)
        time<-as.numeric(colnames(surv))
        rownames(surv)<-pull(data,.data[[id.var]])
        colnames(surv)<-NULL
        pred_surv(time=time,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                model<-do.call(survivalmodels::coxtime,arg)
                prediction<-predict(model,data%>%filter(.data[[id.var]] %in% .env$fold),type="survival",distr6=FALSE)
                time<-as.numeric(colnames(prediction))
                colnames(prediction)<-NULL
                surv<-lapply(all.times,function(t){
                    i<-find.first.TRUE.index(time<=t,noTRUE=0)
                    if(i==0){
                        out<-matrix(1,nrow=length(fold),ncol=1)
                    }else{
                        out<-prediction[,i]
                    }
                    out
                })%>%do.call(what=cbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}





#' @title Wrapper of `survivalmodels::deepsurv`
#' @name fit_deepsurv
#' @param formula formula used by \code{\link[survivalmodels:deepsurv]{survivalmodels::deepsurv}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survivalmodels:deepsurv]{survivalmodels::deepsurv}}. We encourage using a named list. Will be passed to \code{\link[survivalmodels:deepsurv]{survivalmodels::deepsurv}} by running a command like `do.call(deepsurv, option)`. The user should not specify `formula`, `data` and `reverse`; `time_variable`, `status_variable`, `x`, `y` will be ignored.
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_deepsurv<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),...){
    .requireNamespace("survivalmodels")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data","reverse") %in% names(option))){
        stop("option specifies formula, data or reverse")
    }
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        model<-do.call(survivalmodels::deepsurv,arg)
        surv<-predict(model,type="survival",distr6=FALSE)
        time<-as.numeric(colnames(surv))
        rownames(surv)<-pull(data,.data[[id.var]])
        colnames(surv)<-NULL
        pred_surv(time=time,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                model<-do.call(survivalmodels::deepsurv,arg)
                prediction<-predict(model,data%>%filter(.data[[id.var]] %in% .env$fold),type="survival",distr6=FALSE)
                time<-as.numeric(colnames(prediction))
                colnames(prediction)<-NULL
                surv<-lapply(all.times,function(t){
                    i<-find.first.TRUE.index(time<=t,noTRUE=0)
                    if(i==0){
                        out<-matrix(1,nrow=length(fold),ncol=1)
                    }else{
                        out<-prediction[,i]
                    }
                    out
                })%>%do.call(what=cbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}





#' @title Wrapper of `survivalmodels::dnnsurv`
#' @name fit_dnnsurv
#' @param formula formula used by \code{\link[survivalmodels:dnnsurv]{survivalmodels::dnnsurv}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survivalmodels:dnnsurv]{survivalmodels::dnnsurv}}. We encourage using a named list. Will be passed to \code{\link[survivalmodels:dnnsurv]{survivalmodels::dnnsurv}} by running a command like `do.call(dnnsurv, option)`. The user should not specify `formula`, `data` and `reverse`; `time_variable`, `status_variable`, `x`, `y` will be ignored.
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_dnnsurv<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),...){
    .requireNamespace("survivalmodels")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data","reverse") %in% names(option))){
        stop("option specifies formula, data or reverse")
    }
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        model<-do.call(survivalmodels::dnnsurv,arg)
        surv<-predict(model,type="survival",distr6=FALSE)
        time<-as.numeric(colnames(surv))
        rownames(surv)<-pull(data,.data[[id.var]])
        colnames(surv)<-NULL
        pred_surv(time=time,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                model<-do.call(survivalmodels::dnnsurv,arg)
                prediction<-predict(model,data%>%filter(.data[[id.var]] %in% .env$fold),type="survival",distr6=FALSE)
                time<-as.numeric(colnames(prediction))
                colnames(prediction)<-NULL
                surv<-lapply(all.times,function(t){
                    i<-find.first.TRUE.index(time<=t,noTRUE=0)
                    if(i==0){
                        out<-matrix(1,nrow=length(fold),ncol=1)
                    }else{
                        out<-prediction[,i]
                    }
                    out
                })%>%do.call(what=cbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}




#' @title Wrapper of `survivalmodels::akritas`
#' @name fit_akritas
#' @param formula formula used by \code{\link[survivalmodels:akritas]{survivalmodels::akritas}}
#' @param data data containing all covariates, follow-up time, event indicator and id
#' @param id.var see \code{\link{SDRsurv}}
#' @param time.var see \code{\link{SDRsurv}}
#' @param event.var see \code{\link{SDRsurv}}
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to \code{\link[survivalmodels:deepsurv]{survivalmodels::akritas}}. We encourage using a named list. Will be passed to \code{\link[survivalmodels:akritas]{survivalmodels::akritas}} by running a command like `do.call(akritas, option)`. The user should not specify `formula`, `data` and `reverse`; `time_variable`, `status_variable`, `x`, `y` will be ignored.
#' @param lambda bandwidth parameter for uniform smoothing kernel in nearest neighbours estimation. The default value of 0.5 is arbitrary and should be chosen by the user
#' @param ... ignored
#' @return a \code{\link{pred_surv}} class containing fitted survival curves for individuals in `data`
#' @export
fit_akritas<-function(formula,data,id.var,time.var,event.var,nfold=1,option=list(),lambda=0.5,...){
    .requireNamespace("survivalmodels")
    
    #check if option is a list and whether it specifies formula and data
    assert_that(is.list(option))
    if(any(c("formula","data","reverse") %in% names(option))){
        stop("option specifies formula, data or reverse")
    }
    
    if(all(pull(data,.data[[event.var]])==0)){
        return(fit_no_event(data,id.var))
    }
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        arg<-c(
            list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
            option
        )
        model<-do.call(survivalmodels::akritas,arg)
        surv<-predict(model,type="survival",lambda=lambda,distr6=FALSE)
        time<-as.numeric(colnames(surv))
        rownames(surv)<-pull(data,.data[[id.var]])
        colnames(surv)<-NULL
        pred_surv(time=time,surv=surv)
    }else{
        folds<-create.folds(pull(data,.data[[id.var]]),nfold)
        surv.list<-lapply(folds,function(fold){
            if(data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%pull(.data[[event.var]])%>%{all(.==0)}){
                surv<-matrix(1,nrow=length(fold),ncol=length(all.times))
                rownames(surv)<-fold
                surv
            }else{
                arg<-c(
                    list(formula=formula,data=data%>%filter(!(.data[[id.var]] %in% .env$fold))%>%select(!.data[[id.var]])), #remove id.var to allow for . in formula
                    option
                )
                model<-do.call(survivalmodels::akritas,arg)
                prediction<-predict(model,data%>%filter(.data[[id.var]] %in% .env$fold),type="survival",lambda=lambda,distr6=FALSE)
                time<-as.numeric(colnames(prediction))
                colnames(prediction)<-NULL
                surv<-lapply(all.times,function(t){
                    i<-find.first.TRUE.index(time<=t,noTRUE=0)
                    if(i==0){
                        out<-matrix(1,nrow=length(fold),ncol=1)
                    }else{
                        out<-prediction[,i]
                    }
                    out
                })%>%do.call(what=cbind)
                rownames(surv)<-fold
                surv
            }
        })
        surv<-do.call(rbind,surv.list)
        surv<-surv[order(rownames(surv)),,drop=FALSE]
        pred_surv(time=all.times,surv=surv)
    }
}
