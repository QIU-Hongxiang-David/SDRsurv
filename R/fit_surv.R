#' @title Options of machine learning methods' wrappers for fitting conditional survival curves
#' @name fit_surv_option
#' @param nfold number of folds used when fitting survival curves with sample splitting. Default is 1, meaning no sample splitting
#' @param option a list containing optional arguments passed to the wrapped machine learning function. Will be used in a command like `do.call(machine.learning, option)`. `formula` and `data` should not be specified. For \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}, if `tune=TRUE`, then `mtry` and `nodesize` should not be specified either.
#' @param oob whether to use out-of-bag (OOB) fitted values from random forests (\code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} and \code{\link[party:cforest]{party::cforest}}) when sample splitting is not used (`nfold=1`). Ignored otherwise.
#' @param tune whether to tune `mtry` and `nodesize` for \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}}. Ignored for other methods.
#' @param tune.option a list containing optional arguments passed to \code{\link[randomForestSRC:tune]{randomForestSRC::tune.rfsrc}} if \code{\link[randomForestSRC:rfsrc]{randomForestSRC::rfsrc}} is used and `tune=TRUE`; ignored otherwise. `doBest` should not be specified.
#' @export
fit_surv_option<-function(nfold=1,option=list(),oob=TRUE,tune=TRUE,tune.option=list()){
    assert_that(is.count(nfold))
    assert_that(is.flag(oob))
    assert_that(is.flag(tune))
    out<-list(nfold=nfold,option=option,oob=oob,tune=tune,tune.option=tune.option)
    class(out)<-"fit_surv_option"
    out
}

fit_surv<-function(method=c("rfsrc","ctree","rpart","cforest","coxph"),...){
    method<-match.arg(method)
    if(method=="rfsrc"){
        fit_rfsrc(...)
    }else if(method=="ctree"){
        fit_ctree(...)
    }else if(method=="rpart"){
        fit_rpart(...)
    }else if(method=="cforest"){
        fit_cforest(...)
    }else if(method=="coxph"){
        fit_coxph(...)
    }
}


fit_no_event<-function(data,id.var){
    surv<-matrix(1,nrow=nrow(data),ncol=1)
    rownames(surv)<-pull(data,.data[[id.var]])
    pred_surv(Inf,surv)
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
    
    if(nfold==1){
        if(all(pull(data,.data[[event.var]])==0)){
            fit_no_event(data,id.var)
        }else{
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
        }
    }else{
        all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
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
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        if(all(pull(data,.data[[event.var]])==0)){
            fit_no_event(data,id.var)
        }else{
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
        }
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
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        if(all(pull(data,.data[[event.var]])==0)){
            fit_no_event(data,id.var)
        }else{
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
        }
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
    
    all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
    if(nfold==1){
        if(all(pull(data,.data[[event.var]])==0)){
            fit_no_event(data,id.var)
        }else{
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
        }
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
    if(nfold==1){
        if(all(pull(data,.data[[event.var]])==0)){
            fit_no_event(data,id.var)
        }else{
            arg<-c(
                list(formula=formula,data=select(data,!.data[[id.var]])), #remove id.var to allow for . in formula
                option
            )
            cox.model<-do.call(survival::coxph,arg)
            model<-survival::survfit(cox.model,newdata=data)
            surv<-t(as.matrix.rowvec(model$surv))
            rownames(surv)<-pull(data,.data[[id.var]])
            pred_surv(time=model$time,surv=surv)
        }
    }else{
        all.times<-data%>%filter(.data[[event.var]]==1)%>%pull(.data[[time.var]])%>%unique%>%sort
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
                model<-survival::survfit(cox.model,newdata=data%>%filter(.data[[id.var]] %in% .env$fold))
                model.surv<-t(as.matrix.rowvec(model$surv))
                surv<-lapply(all.times,function(t){
                    i<-find.first.TRUE.index(model$time.interest<=t,noTRUE=0)
                    if(i==0){
                        out<-matrix(1,nrow=length(fold),ncol=1)
                    }else{
                        out<-as.matrix(model$surv)[,i]
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
