#' @title Prepare data used in regression for estimating (conditional) suvival probabilities
#' @name prep.Q.reg.data
#' @description Sum up pseudo-outcomes in each time window and produce the final pseudo-outcome in the regression to estimate P(T <= t | T > truncation time, covariates available at truncation time); link covariates at each check-in time up to truncation index to form the data with all historic covariates that can be used in the regression. The output can be used by \code{\link[SuperLearner]{SuperLearner}} or other regression methods for the final estimation.
#'
#' @param covariates See \code{\link{SDRsurv}}.
#' @param stagewise.pseudo.outcomes a list of list of data frames. Each element in the outer layer corresponds to a value t in `tvals` for which to estimate the conditional survival probability. Each element in the inner layer corresponds to a data frame at a check-in time (i.e., stage) before t. Each data frame should contain a variable to identify each individual (see `id.var` argument of \code{\link{SDRsurv}}) and a variable named `pseudo.outcome` that contains the pseudo-outcome (e.g., doubly robust transformation) for that check-in time. For check-in times before the truncation time, `pseudo.outcome` should equal to 0.
#' @param truncation.index index of the check-in time to which left-truncation is applied. The truncation time is `check.in.times[truncation.index]`. Covariates available up to (inclusive) `check.in.times[truncation.index]` are of interest. See \code{\link{SDRsurv}}
#' @param id.var name of the variable that identifies each individual. Just provide the expression such as \code{id} and do not embrace it with `""` or `''`.
#' @return A list containing:
#' \enumerate{
#' \item `pseudo.outcomes`: a list of data frames with pseudo-outcomes. Each data frame corresponds to one t in `tvals` (see \code{\link{SDRsurv}}) and contains two variables: a variable with same name as `id.var`; a variable named `pseudo.outcome`
#' \item `history`: a data frame containing all historic variables
#' }
#' All data frames are sorted by `id.var`.
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{SDRsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @export
prep.Q.reg.data<-function(
    covariates,
    stagewise.pseudo.outcomes,
    truncation.index,
    id.var
){
    all.have.pseudo.outcome<-all(sapply(stagewise.pseudo.outcomes,function(Y.list){
        all(sapply(Y.list,has_name,which="pseudo.outcome"))
    }))
    if(!all.have.pseudo.outcome){
        stop("1+ data frame in stagewise.pseudo.outcomes does not have variable pseudo.outcome")
    }
    
    all.have.id.var<-all(sapply(stagewise.pseudo.outcomes,function(Y.list){
        all(sapply(Y.list,has_name,which=id.var))
    }))
    if(!all.have.id.var){
        stop(paste("1+ data frame in stagewise.pseudo.outcomes does not have",id.var))
    }
    
    
    history<-reduce(covariates[1:truncation.index],.f=function(d1,d2){
        right_join(d1,d2,by=id.var)
    })%>%arrange(.data[[id.var]])

    pseudo.outcomes<-lapply(stagewise.pseudo.outcomes,function(Y.list){
        reduce(Y.list,.f=function(d1,d2){
            d2<-d2%>%rename(pseudo.outcome2=.data$pseudo.outcome)
            left_join(d1,d2,by=id.var)%>%
                mutate(pseudo.outcome2=ifelse(is.na(.data$pseudo.outcome2),0,.data$pseudo.outcome2),
                       pseudo.outcome=.data$pseudo.outcome+.data$pseudo.outcome2)%>%
                select(!.data$pseudo.outcome2)
        })%>%arrange(.data[[id.var]])
    })
    
    list(pseudo.outcomes=pseudo.outcomes,history=history)
}




#' @title Estimate (conditional) survival probabilities given pseudo-outcomes with SuperLearner
#' @name estQ.SuperLearner
#' @description Estimate P(T <= t | T > truncation time, covariates available at truncation time) for given t, where T is the time to event.
#'
#' @param covariates See \code{\link{SDRsurv}}.
#' @param stagewise.pseudo.outcomes a list of list of data frames. Each element in the outer layer corresponds to a value t in `tvals` for which to estimate the conditional survival probability. Each element in the inner layer corresponds to a data frame at a check-in time before t. Each data frame should contain a variable to identify each individual (see `id.var` argument of \code{\link{SDRsurv}}) and a variable named `pseudo.outcome` that contains the pseudo-outcome (e.g., doubly robust transformation) for that check-in time. For check-in times before the truncation time, `pseudo.outcome` should equal to 0.
#' @param truncation.index index of the check-in time to which left-truncation is applied. The truncation time is `check.in.times[truncation.index]`. Covariates available up to (inclusive) `check.in.times[truncation.index]` are of interest. See \code{\link{SDRsurv}}
#' @param id.var (character) name of the variable that identifies each individual.
#' @param Q.formula formula to specify covariates being used for estimating P(T <= t | T > `check.in.times[truncation.index]`, covariates available at `check.in.times[truncation.index]`). Set to include intercept only (`~ 0` or `~ -1`) for marginal survival probability, which is simply the mean of pseudo-outcomes. Default is `~ .`, which includes main effects of all available covariates up to (inclusive) the `truncation.time`.
#' @param Q.SuperLearner.control a list containing optional arguments passed to \code{\link[SuperLearner]{SuperLearner}}. We encourage using a named list. Will be passed to \code{\link[SuperLearner]{SuperLearner}} by running a command like `do.call(SuperLearner, Q.SuperLearner.control)`. The user should not specify `Y`, `X` and `family`. This argument is ignored if `Q.formula` has no covariates.
#' @return a list of fitted `SuperLearner` models corresponding to each t in `tvals`. If `Q.formula` is empty, then return a list of numbers, each being estimated P(T <= t) for t in `tvals`.
#' @section Warning:
#' This function is designed to be called by other functions such as \code{\link{SDRsurv}}, therefore inputs are not thoroughly checked. Incorrect inputs may lead to errors with non-informative messages. The user may call this function if more flexibility is desired.
#' @section Custom learners:
#' Custom learners may be specified by providing an element named `SL.library` in `Q.SuperLearner.control`.The user may refer to resources such as \url{https://cran.r-project.org/web/packages/SuperLearner/vignettes/Guide-to-SuperLearner.html} for a guide to create custom learners.
#' @export
estQ.SuperLearner<-function(
    covariates,
    stagewise.pseudo.outcomes,
    truncation.index,
    id.var,
    Q.formula=~.,
    Q.SuperLearner.control=list(SL.library="SL.lm")
){
    prepared<-prep.Q.reg.data(covariates,stagewise.pseudo.outcomes,truncation.index,id.var)
    pseudo.outcomes<-prepared$pseudo.outcomes
    history<-prepared$history
    
    assert_that(is.list(Q.SuperLearner.control))
    if(any(c("Y","X","family") %in% names(Q.SuperLearner.control))){
        stop("Q.SuperLearner.control specifies Y, X or family!")
    }
    
    X<-model.frame(Q.formula,data=history%>%arrange(.data[[id.var]])%>%select(!.data[[id.var]]))
    lapply(pseudo.outcomes,function(d){
        d<-arrange(d,.data[[id.var]])
        Y<-d$pseudo.outcome
        if(ncol(X)==0){ #empty model
            mean(Y)
        }else{
            SuperLearner.arg<-c(
                list(Y=Y,X=X,family=gaussian()),
                Q.SuperLearner.control
            )
            do.call(SuperLearner,SuperLearner.arg)
        }
    })
}
