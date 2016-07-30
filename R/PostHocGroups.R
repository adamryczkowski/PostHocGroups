#' Provides algorithm that divides ordered sets of items into easy to interpret exclusive groups based on
#' binary relation matrix which tells whether one item is different from another or more general
#' dissimilarity matrix, which also tells how significantly one item is different from another.
#'
#'
#' One usefull application is to present results of pairwise post-hoc tests of many groups in a
#' easy to understand manner.
#'
#' This matrix is input either explicitely as \code{dissimilarity.matrix} or inferred implicitely
#' from input post-hoc model.
#'
#' A given dissimilarity matrix is not divisable into the exclusive groups unless we accept some error. The algorithm
#' operates by trying to maximize quality of partitioning (QoP).
#'
#' The algorithm doesn't search through all possible ways to divide the item set -
#' the number of iterations would equal \eqn{2^n} (and this value do not include recursively divided groups ) and
#' even for relatively small number of post-hoc groups it would be computationaly not feasible to calculate.
#' It rather uses hierarchical clustering to do the heavy lifting. Details in the \strong{algorithm} section.
#'
#' @section Features:
#'
#'\itemize{
#'
#' \item The algorithm can output alternative results for user to choose among.
#'
#' \item Algorithm can identify a case, when certain group of items is not different from all (almost all)
#' other exclusive groups and output it separately.
#' This can happen when performing post hoc tests on grouping variable that doesn't divide the population equally.
#'
#' \item Algorithm support nested exclusive groups, i.e. it can identify a case, when there is a subgroup of variables
#' that are together different from all other groups, but are not all similar to each other. Nested grouping
#' can only improve fit to the given dissimilarity matrix.
#'}
#'
#'
#'
#'
#'
#'
#' @section The algorithm:
#'
#' Basically the algorithm tries to get as many as possible meaningful partitionings of the items and
#' chooses the one which maximises the utility.
#'
#' \subsection{Utility of the partitioning}{
#'
#' The partition to evalute is converted into the Result Dissimilarity Matrix (RDM), where situation where two items
#' belong to the same partition is marked as "0" and if the two items belong to different partitions - by "1".
#' The utility of the paritioning is simply a difference between sum of all elements of the original
#' dissimilarity matrix where the corresponding elements of the RDM matrix have value "0" minus the elements, where
#' corresponding elements of the RDM matrix have value "1".
#'
#' In the following description for simplicity I will assume that the dissimilarity matrix is binary, i.e. it tells
#' whether two items are different from each other (value 1) as opposed to being similar (value -1).
#'
#' I will also asume, that the items are sorted according to the same measure that user used to calculate
#' the dissimilarity matrix (usually it is group mean in case of PostHoc tests).
#'
#'}
#'
#' \subsection{The algorithm step-by-step}{
#'
#' The algorithm tries very hard not to test all \eqn{2^n} possible partitionings.
#'
#' \enumerate{
#'
#' \item Pre-processing of the dissimilarity matrix
#'
#' At the beginning the algorithm will check if user input binary dissimilarity matrix
#' (i.e. if it contains only 2 distinct values: is equal/is nonequal). If so, then a random noise will be added.
#' The noise will be sufficiently small to not to change the meaning of the matrix. It is needed because a lot
#' of similar values in the dissimilarity matrix would lead to sub-optimal solutions.
#'
#' \item Gathering all possible partitionings of the items
#'
#' \enumerate{
#' \item \strong{Getting the "smile vector"}:
#' The key statistics to get meaningful candidates for a common group is the "smile vector". It is a
#' column-wise sum of the dissimilarity matrix (or row-wise: it doesn't matter, since the
#' matrix is assumed symmetric) which can be interpreted as a number of significant differences for each item.
#' Elements of the "smile vector" correspond to items that are going to partitioned.
#' Usually the edge items (those with minimal of maximal value)
#' tend to be more dissimilar to the rest of items, than items that have value closer to the median.
#' If the vector is plotted it usually has a "U" shape, hence its name: the "smile vector".
#'
#' \item \strong{Correcting the "smile vector" for its smile}:
#' [This is optional step. It is performed only on top recursion level, and only if user ask for it]
#' To offset the "U" shape of the "smile vector", algorithm can fits second level multiplicative
#' "penalty" polynomial in the form
#' \deqn{p(x) = A ( \frac{2 x^2}{(n-1)^2}-\frac{2 x}{n-1}+1 )}.
#' The polynomial is choosen the way, that its second derivative is always positive, and \eqn{p(0)=A},
#' \eqn{p(n-1)=A} and \eqn{p(\frac{n-1}{2})=\frac{A}{2}}. In other words, the penalty makes sure, that significant
#' differences between items that are on the far ends of the item list are twice less pronounced than differences
#' in the middle of the set.
#' After fitting the only parameter of the penalty (\eqn{A}) it multiplies all elements of the smile vector by the
#' fitted function to get the corrected smile vector.
#'
#' \item Partition the smile vector into the common group and the items that will be further partitioned.
#'
#' Common group is a group of items that are not dissimilar to any other items in the set. Items in this group
#' would have very low sum of dissimilarities to the other items. Algorithm guesses the common group by assigning
#' the first \eqn{i \in (0,n)} elements of the sorted (ascending) smile vector.
#'
#' \item Perform cluster analysis on the items that are left
#'
#' Partitioning is performed by the hierarchical clustering \code{hclust} using the dissimilarity matrix direclty.
#'
#' Algorithm gathers for evaluation all partitionings for each distinct "height" of the hierarchical clustering.
#'
#' }
#' \item If the final recursion level is reached, return all alternatives together with the utility (quality)
#'
#' \item Otherwise, prepare for efficient nested partitioning of all eligible groups.
#'
#'  The algorithm is complicated, because of the need to restate the problem in manner compatible with the main
#'  function, and to make sure we don't launch the recursion for the same set of items twice.
#' }
#' }
#' @docType package
#'
#' @name PostHocGroups
NULL

# nocov start
.onLoad	<-	function(libname,	pkgname)	{
  invisible()
}
# nocov end

#' Groups result from post-hoc test
#'
#' The function is generic and can consume either dissimilarity matrix (of class \code{dist},
#' e.g. produced by \code{stats::dist} function), or
#' \code{glht} object, produced by the \code{\link[multcomp]{glht}} function
#' from the \code{multcomp} package.
#'
#'
#' @param object Object of one of the supported types.
#' @param max.recursive.level Defaults to 1 for flat representation.
#'        If more, algorithm will try to recursively subpartition groups that have
#'        more than 2 elements.
#' @param solutions.count Number of sultions to report. Defaults to 1, i.e. only the best solution
#'        will be printed.
#' @return Object of class \code{PostHocGroups}. If only one possible solution was requested,
#' the object has also class \code{PostHocGroup}.
#' @export
GroupPostHocs<-function(object,...)
{
  UseMethod('GroupPostHocs')
}

#' Groups result from the \code{glht} post-hoc test
#'
#' This is a specialization of the generic function \code{GroupPostHocs}
#'
#' @param res Object of \code{glht} class.
#' @param max.recursive.level Defaults to 1 for flat representation.
#'        If more, algorithm will try to recursively subpartition groups that have
#'        more than 2 elements.
#' @param solutions.count Number of sultions to report. Defaults to 1, i.e. only the best solution
#'        will be printed.
#' @return Object of class \code{PostHocGroups}. If only one possible solution was requested,
#' the object has also class \code{PostHocGroup}.
#' @export
GroupPostHocs.glht<-function(res, max.recursive.level=1, solutions.count=1)
{
  gr<-rownames(res$model$contrasts$varnr)
  symmetric.list.outer2<-function(res)
  {
    len<-length(gr)
    pvalues<-res$test$pvalues

    pairs<-plyr::aaply(.data=names(res$test$tstat), .fun = function(x, gr)
    {
      x2<-stringr::str_match_all(x, '^([^\\s^\\-]+)\\s-\\s([^\\s^\\-]+)$')[[1]][-1]
      return(which(gr %in% x2))
    }, .margins=1, gr=gr
    )

    out<-matrix(nrow=len,ncol=len)
    colnames(out)<-gr
    rownames(out)<-gr
    for (k in seq_along(pvalues))
    {
      i<-pairs[[k,1]]
      j<-pairs[[k,2]]
      tmp<- min(-log(pvalues[[k]])+log(0.05),20)
      out[[i,j]]<-tmp
      out[[j,i]]<-tmp
    }
    return(out)
  }
  means<-plyr::laply(gr, function(x){mean(res$model$data[varnr==x,value], na.rm=TRUE)})
  m<-symmetric.list.outer2(res)
  row<-GroupPostHocs.dist(dissimiliarity.matrix = m,
                          means = means,
                          max.recursive.level = max.recursive.level)
  return(row)
}

#' Groups elements of the given \code{dissimilarity matrix}
#'
#' This is the most general way to invoke the algorithm. User needs to provide distance matrix
#' that is interpreted as dissimilarity matrix, with value that is interpreted as a
#' measure of dissimilarity between items:
#' \describe{
#'  \item{matrix value >0}{Two items are similar}
#'  \item{matrix value <0}{Two items are different}
#' }
#'
#' To get meaningfull results, user should either provide matrix with items (columns/rows)
#' of the matrix are sorted according to the corresponding value of the tested item, or
#' provide the \code{means} argument and the algorithm will sort the elements internally.
#'
#' @param dissimiliarity.matrix Object of \code{dist} class.
#' @param means Optional parameter. Requires numeric vector of the size of number of columns of
#'        the dissimiliarity.matrix or \code{NULL} (default). If specified, function will sort
#'        the dissimiliarity.matrix according to those values during the pre-process stage of the
#'        algorithm.
#' @param max.recursive.level Defaults to 1 for flat representation.
#'        If more, algorithm will try to recursively subpartition groups that have
#'        more than 2 elements.
#' @param solutions.count Number of sultions to report. Defaults to 1, i.e. only the best solution
#'        will be printed.
#' @return Object of class \code{PostHocGroups}. If only one possible solution was requested,
#' the object has also class \code{PostHocGroup}.
#' @export
GroupPostHocs.dist<-function(dissimiliarity.matrix, means=NULL, max.recursive.level=1, solutions.count=1)
{
  myprox<-as.matrix(dissimiliarity.matrix)
  mylen<-ncol(myprox)
  solutions<-data.table::as.data.table(read.table(text="",colClasses = 'character',col.names = 1:mylen))
  colnames(solutions)<-as.character(c(1:mylen))

  solutions<-recursive.clustering(myprox=myprox,
                       prefix='',
                       solutions=solutions,
                       recursive.level=1,
                       max.recursive.level=max.recursive.level,
                       means=means,
                       flag.smile.correct=TRUE)

  rows<-solutions[order(-qual)][1:solutions.count,!colnames(solutions) %in%
                                    c('id','sum.inner','sum.outer','n.inner','n.outer'),with=FALSE]
  #browser()
  rows2<-matrix(as.character(rows),nrow=nrow(rows),ncol=ncol(rows))
  colnames(rows2)<-c(colnames(myprox),c('_qual'))
  class(rows2)<-'PostHocGroups'

  if (solutions.count==1)
    class(rows2)<-c('PostHocGroup',class(rows2))
  return(rows2)
}

correct.for.smile<-function(smile)
{
  #browser()
  mylen=length(smile)
  f<-function(x)
  { #Function with "ideal" smile: f(0)=1, f(mylen-1)=1, f((mylen-1)/2)=0.5
    return((2*x*x)/((mylen-1)*(mylen-1))-(2*x)/(mylen-1)+1)
  }
  x<-seq(0,mylen-1)

  A <- tryCatch({
    #fitted the only degree of freedom: amplitude of the function.
    coef(nls(smile ~ A*f(x), start=list(A=mean(smile))))
  }, error = function(e) {
    1
  })
  smile2<-plyr::aaply(x,1,function(x){f(x)*A}) #Idealistic, fitted smile
  return(smile-smile2) #Corrected smile.
}

add.possible.partitionings<-function(myprox, solutions, solution.idx, prefix)
{
  #browser()
  htree<-hclust(d=as.dist(myprox))
  heights<-unique(htree$height) #These heights enumerate all possible partitions of the items into exclusive groups
  row<-rep(paste0(prefix, '0'), ncol(solutions))
  plyr::a_ply(heights,1,
              function(height)
              {
                t<-cutree(tree=htree,k=NULL,h=height)
                #browser()
                if (length(row)>length(solution.idx) && #We consider common group items and...
                    length(unique(t))<=1) #...and there is only on exclusive group...
                {
                  return()  #...don't add this rubbish to the solutions database.
                }
                row[solution.idx]<-paste0(prefix,t)
                solutions<<-append.row(solutions, row)
              })
  row<-rep(paste0(prefix, '0'), ncol(solutions))
  row[solution.idx]<-paste0(prefix, seq_along(solution.idx))
  solutions<-append.row(solutions, row) #Every element belongs to distinct group
  return(solutions)
}

#' This is the main function.
recursive.clustering<-function(myprox,
                               prefix,
                               solutions,
                               recursive.level,
                               max.recursive.level,
                               means=NULL,
                               flag.smile.correct=TRUE,
                               flag.find.common.set=TRUE)
{
  if (length(unique(as.numeric(myprox)))<5)
  {
    #browser()
    #we have binary dissimilarity matrix. Need to add some tiny randomness to it
    threshold<-1/(ncol(myprox)*4) #That is about as much as we can. If we add more, we risk significant change
    #to the dissimilarity matrix.
    salt<-matrix(data=0, ncol=ncol(myprox),nrow=nrow(myprox))
    salt[lower.tri(salt)]<-runif(sum(lower.tri(myprox)), min=-threshold, max=threshold)
    salt<-salt+t(salt)
    myprox.salted<-myprox+salt #myprox.salted is used for clustering only, while original myprox is used
    # for calculating quality
  } else {
    myprox.salted<-myprox
  }
  if (is.null(colnames(myprox)))
  {
    colnames(myprox)<-as.character(seq(1,ncol(myprox)))
  }
#  browser()

  if (flag.find.common.set)
  {
    sig.count<-plyr::aaply(myprox.salted,2,sum, na.rm=TRUE) #vector with total number of significant differences
    # for each item. We will do 1-dim clustering on this vector to separate common group from exclusive groups.
    mylen<-ncol(myprox)
    if (flag.smile.correct)
    {
      if (!is.null(means))
      {
        ord<-order(means)
        #Items in the middle should have half the
        #amount of significant differences, than variables on the edges, that's why the name.
        smile<-sig.count[ord]
      } else {
        smile<-sig.count
        ord<-seq(1,mylen)
      }
      sig.count.cor<-correct.for.smile(smile)[Matrix::invPerm(ord)]
    } else {
      sig.count.cor<-sig.count
    }
    commonSetPartitions<-sort(unique(sig.count.cor))
  } else {
    commonSetPartitions<- -Inf
  }

  for(i in seq_along(commonSetPartitions)) #Number of possible divisions into clusters into which we divide all items. We need theese clusters
    #to get many meaningful divisions of the set of into 2 groups: one that gets divided into similar clusters
    #based on user-supplied dissimilarity matrix, and one that is ignored.
  {
    h<-commonSetPartitions[[i]] #
#    common.items<-colnames(myprox)[sig.count.cor<h] #Those items are left behind because
    # have small number of dissimilarities to other items
    partition.items<-colnames(myprox)[sig.count.cor>=h]
    # those items later will be partitioned into exclusive sets.
    if (length(partition.items)>1)
    {
      idx<-which(colnames(myprox) %in% partition.items)
#      myprox2<-as.dist(myprox.salted[idx,idx])

      solutions<-add.possible.partitionings(myprox = myprox.salted[idx,idx],
                                           solutions = solutions,
                                           solution.idx = idx,
                                           prefix = prefix)

    }
  }

  # Stripping not used records from solutions
  n<-attr(solutions, 'rowcount')
  #browser()
  solutions<-data.table::as.data.table(solutions[1:n,])
  # Removing duplicate records
  solutions<-unique(solutions)

  if (recursive.level < max.recursive.level)
  {
    browser()
    #We need to restate problem into a from that
    # a) is recursable in the current function
    # b) doesn't launch duplicate recrsions
    podzbiory.hash<-hash::hash() #This hash is a substitute for memoization of the recursive calls
    sol2podzbiory<-list() #Solution ID -> list(podzbior.id,...)
    j<-1
    for (i in 1:nrow(solutions))
    {
      row<-as.character(solutions[i,])
      clusters<-table(row)
      clusters<-clusters[names(clusters)!='0' & clusters>2] #We remove all groups that consist of less than 2 items.
      #These groups are not eligible for recursive checks, because there is nothing to gain
      podzbiory.list<-list()
      plyr::a_ply(names(clusters),1,function(cl)
        {
          r<-which(row %in% cl)
          h<-paste0(r,collapse=',')
          podzbiory.list<<-c(podzbiory.list, list(list(class=cl, hash=h)))
          podzbiory.hash[[h]]<-TRUE
          j<<-j+1
        })
      sol2podzbiory[[i]]<-podzbiory.list
    }
    browser()
    podzbiory<-hash::keys(podzbiory.hash)
    podsolutions.set<-list()
    for(i in seq_along(podzbiory))
    {
      podzbior<-as.integer(unlist(strsplit(podzbiory[[i]],',')))
      podmyprox<-myprox[podzbior,podzbior]
      colnames(podmyprox)<-podzbior
      rownames(podmyprox)<-podzbior
      podmylen<-length(podzbior)
      podsolutions<-data.table::as.data.table(read.table(text="",colClasses = 'character',col.names = 1:podmylen))
      colnames(podsolutions)<-as.character(podzbior)
      podmeans<-means[podzbior]

      podsolutions<-recursive.clustering(myprox = podmyprox,
                                         prefix = '',
                                         solutions = podsolutions,
                                         recursive.level = recursive.level+1,
                                         max.recursive.level = max.recursive.level,
                                         means = podmeans,
                                         flag.smile.correct = FALSE,
                                         flag.find.common.set = flag.find.common.set)
      if (!is.null(podsolutions) & nrow(podsolutions)>0)
      {
        podsolutions.set[[as.character(i)]]<-podsolutions
      }
    }

    browser()
    solutions[,id:=seq(1,nrow(solutions))]
    solutions[,c("sum.inner","sum.outer", "n.inner", "n.outer", "qual"):=-1]
    qualCols<-which(colnames(solutions) %in% c('sum.inner','sum.outer','n.inner','n.outer','qual'))

    for (i in as.integer(1:nrow(solutions)))
    {
      solpodzbiory<-sol2podzbiory[[i]]
      complex.podzbiory<-integer(0)
      complex.classes<-character(0)
      for (j in seq_along(solpodzbiory))
      {
        solpodzbior<-solpodzbiory[[j]]
        h<-solpodzbior$hash
        hid<-which(podzbiory==h)
        podsolutions<-podsolutions.set[[as.character(hid)]]
        if (!is.null(podsolutions))
        {
          quals<-podsolutions[['qual']]
          bestqual<-which(quals==max(quals))
          n.outer<-podsolutions[bestqual,n.outer]
          if (n.outer>0)
          {
            complex.podzbiory<-c(complex.podzbiory,j)
            complex.classes<-c(complex.classes, solpodzbior$class)
          }
        }
      }

      cat(paste0('i = ',i,'\n'))
      #browser()
      solution<-solutions[i,seq(1,ncol(myprox)),with=FALSE]
      qual<-qual.of.partition(
        classification.matrix = as.matrix(makeprox(as.character(solution))),
        dissimilarity.matrix = myprox,
        class.to.ignore = complex.classes)

      data.table::set(solutions, i, qualCols[[1]], qual$inner.util)
      data.table::set(solutions, i, qualCols[[2]], qual$outer.util)
      data.table::set(solutions, i, qualCols[[3]], qual$n.inner.util)
      data.table::set(solutions, i, qualCols[[4]], qual$n.outer.util)
      data.table::set(solutions, i, qualCols[[5]], qual$qual)
#      solqual<-qual$qual
      if (length(complex.podzbiory)>0)
      {
        #browser()
        cat(paste0(i, ' ma ', length(complex.podzbiory), ' podzbiorów!\n'))
      }

      for (j in complex.podzbiory)
      {
        solpodzbior<-solpodzbiory[[j]]
        h<-solpodzbior$hash
        hid<-which(podzbiory==h)
        podsolutions<-podsolutions.set[[as.character(hid)]]
        quals<-podsolutions[['qual']]
        bestqual<-which(quals==max(quals))
        data.table::set(solutions,i,qualCols[[1]],qual$inner.util + podsolutions[bestqual, sum.inner] + podsolutions[bestqual, sum.outer])
        data.table::set(solutions,i,qualCols[[3]],qual$n.inner.util + podsolutions[bestqual, n.inner] + podsolutions[bestqual, n.outer])
        data.table::set(solutions,i,qualCols[[5]],qual$qual + podsolutions[bestqual, qual])

        solpodzbior<-solpodzbiory[[j]]
        h<-solpodzbior$hash
        podzbior<-as.integer(unlist(strsplit(h,',')))

        newvalues<-as.character(podsolutions[bestqual,as.character(podzbior),with=FALSE])
        solCols<-as.integer(which(colnames(solutions) %in% as.character(podzbior)))
        #browser()
        prefix=paste0(solpodzbiory[[j]]$class,'.')
        for (k in seq_along(solCols))
        {
          data.table::set(solutions, i, solCols[[k]], paste0(prefix,newvalues[[k]]))
        }
      }

    }
    return(solutions)

  } else {
    #browser()
    solutions[,id:=seq(1,nrow(solutions))]
    solutions[,c("sum.inner","sum.outer", "n.inner", "n.outer", "qual"):=
                qual.of.partition(
                  classification.matrix = as.matrix(makeprox(as.character(.SD))),
                  dissimilarity.matrix = myprox),
              by=id,
              .SDcols=seq(1,ncol(myprox))]
  }
  return(solutions)
}

default.util.transform.inner<-function(offset=log(0.05))
{
  return(function(utilities) return(-utilities-offset))
}

default.util.transform.outer<-function(offset=log(0.05))
{
  return(function(utilities) return(utilities+offset))
}

qual.of.partition<-function(classification.matrix, dissimilarity.matrix, class.to.ignore=NULL,
                            util.f.inner.group = default.util.transform.inner(0),
                            util.f.outer.group=default.util.transform.outer(0))
{
  if (is.null(class.to.ignore))
  {
    filtr<-classification.matrix==0 & lower.tri(dissimilarity.matrix)
  } else {
    filtr<-classification.matrix==0 & lower.tri(dissimilarity.matrix)
    #browser()
    filtr[colnames(dissimilarity.matrix) %in% class.to.ignore,colnames(dissimilarity.matrix) %in% class.to.ignore]<-NA
  }
  m0<-util.f.inner.group(dissimilarity.matrix[filtr])
  m0<-sum(c(0,m0), na.rm=TRUE)
  nm0<-sum(filtr,na.rm=TRUE)

  filtr<-classification.matrix==1 & lower.tri(dissimilarity.matrix)
  m1<-util.f.outer.group(dissimilarity.matrix[filtr])
  m1<-sum(c(0,m1), na.rm=TRUE)
  nm1<-sum(filtr, na.rm=TRUE)
  return(list(inner.util=m0, outer.util=m1, n.inner.util=nm0, n.outer.util=nm1, qual=m1+m0))
}

quality.of.fit<-function(row, dissimilarity.matrix)
{
  best.partition<-ifelse(m>0,1,0)
  bestqual<-qual.of.partition(classification.matrix = best.partition,
                              dissimilarity.matrix = dissimilarity.matrix)
  our.partition<-as.matrix(makeprox(row))
  ourqual<-qual.of.partition(classification.matrix = our.partition,
                             dissimilarity.matrix = dissimilarity.matrix)
  qual.difference<-bestqual$qual - ourqual$qual

  diff.partition<-ifelse(best.partition==our.partition,0,ifelse(best.partition-our.partition<0,-1,1))
  type1.errors<-plyr::aaply(diff.partition, 1, function(x){sum(x==1, na.rm=TRUE)})
  type2.errors<-plyr::aaply(diff.partition, 1, function(x){sum(x==-1, na.rm=TRUE)})
  abs.errors<-type1.errors+type2.errors
  sort(type1.errors)
  sort(type2.errors)
  sort(abs.errors)
  return(list(qual.difference=qual.difference,
              type1.errors=type1.errors,
              type2.errors=type2.errors,
              abs.errors=abs.errors))
}
