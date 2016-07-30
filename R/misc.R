#' Creates dissimilarity matrix out from row format got from hierarchical clustering
makeprox<-function(mycutree){
  mylen<-length(mycutree)
  oneitem<-function(x,y)
  {
    if (mycutree[[x]]==0 | mycutree[[y]]==0)
    {
      return (0)
    }
    return ((mycutree[[x]] != mycutree[[y]])*1)
  }
  ans<-array(dim=c(mylen,mylen))
  for(i in 1:mylen)
  {
    for(j in 1:mylen) {
      ans[i,j]<-oneitem(i,j)
    }
  }
  as.dist(ans)
}

#' Returns symmetric simmilarity matrix based on function that returns whether a pair is significant
#'
#' @param items List with items the \code{fun} operates on.
#' @param fn Function with signature \code{fun(item1, item2, ...)} that returns 1 or 0
#'
#' @return Rectangular symmetric matrix with return values from function \code{fn}
#' @export
symmetric.list.outer<-function(items,fun,...)
{
  out<-matrix(nrow=length(items),ncol=length(items))
  for (i in 1:(length(items)-1))
  {
    for (j in (i+1):length(items))
    {
      out[i,j]<-fun(items[[i]],items[[j]],...)
      out[j,i]<-out[i,j]
    }
  }
  out
}

#' Efficiently appends row to the data.table dt.
#' For efficiency comparison of various methods see
#' http://stackoverflow.com/questions/20689650/how-to-append-rows-to-an-r-data-frame/38052208#38052208
append.row<-function(dt, elems)
{
  n<-attr(dt, 'rowcount')
  if (is.null(n))
    n<-nrow(dt)
  if (n==nrow(dt))
  {
    if (n==0)
    {
      m<-4
    } else {
      m<-n
    }
    tmp<-list(rep(NA,m))
    names(tmp)<-colnames(dt)[[1]]
    dt<-data.table::rbindlist(list(dt, tmp) ,fill=TRUE, use.names=TRUE)
    data.table::setattr(dt,'rowcount',n)
  }
  for (j in seq_along(elems))
  {
    data.table::set(dt, i=as.integer(n+1), j, elems[[j]])
  }
  data.table::setattr(dt,'rowcount',n+1)
  return(dt)
}

