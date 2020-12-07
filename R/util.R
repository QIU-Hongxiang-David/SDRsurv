#find the first index of x that is TRUE
#noTRUE: return value when none of x is TRUE
find.first.TRUE.index<-function(x,noTRUE=length(x)+1){
    match(TRUE,x,nomatch=noTRUE)
}

#find the last index of x that is TRUE (loop from the last entry back to the first to and find the first encountered TRUE)
#noTRUE: return value when none of x is TRUE
find.last.TRUE.index<-function(x,noTRUE=0){
    length(x)+1-match(TRUE,rev(x),nomatch=length(x)+1-noTRUE)
}

#' @title Clip all elements in an array to fall in an interval
#' @name clip
#' @param x numeric array
#' @param lower lower end of the interval
#' @param upper upper end of the interval
#' @return clipped array
#' @details For each element `X` in `x`, `clip(X)` equals `lower` if `X`<`lower`, `upper` if `X`>`upper`, and `X` otherwise. This function may be useful to clip predictions to fall in, e.g., \eqn{[0,1]}.
#' @export
clip<-function(x,lower,upper){
    pmax(pmin(x,upper),lower)
}


#' @title Apply administrative censoring to follow up times at a fixed censor.time
#' @name admin.censor
#' @description Right-censor all observations in `follow.up.time` at `censor.time`. This is useful when fitting survival curves within each time window
#' @param follow.up.time See \code{\link{SDRsurv}}
#' @param time.var See \code{\link{SDRsurv}}
#' @param event.var See \code{\link{SDRsurv}}
#' @param censor.time Administrative right-censoring time. Default is `Inf`, i.e., no censoring
#' @return A data frame with the same shape as `follow.up.time` with times and event indicators modified to reflect right-censoring at `censor.time`.
#' @export
admin.censor<-function(follow.up.time,time.var,event.var,censor.time=Inf){
    if(censor.time==Inf){
        follow.up.time
    }else{
        follow.up.time%>%
            mutate(viewed.censored=.data[[time.var]]>.env$censor.time,
                   "{time.var}":=ifelse(.data$viewed.censored,.env$censor.time,.data[[time.var]]),
                   "{event.var}":=ifelse(.data$viewed.censored,0,.data[[event.var]]))%>%
            select(!.data$viewed.censored)
    }
}


#create k folds of a vector id
create.folds<-function(id,k){
    order<-sample.int(length(id))
    d<-suppressWarnings(data.frame(cbind(id[order],1:k)))
    names(d)<-c("id","fold.id")
    lapply(tapply(d$id,d$fold.id,identity,simplify=FALSE),sort)
}
