# Tests a case, when all items belong to one group
cc.all.one.group<-function(item.count=10)
{
  dissimilarity.matrix<-matrix(data = -1, nrow=item.count, ncol=item.count)
  diag(dissimilarity.matrix)<-NA
  pgroups<-PostHocGroups.dist(dissimiliarity.matrix = dissimilarity.matrix)
  best.row<-as.character(pgroups[order(-qual)][1,!colnames(pgroups) %in% c('id','sum.inner','sum.outer','n.inner','n.outer','qual'),with=FALSE])
  testthat::expect_equal(best.row, as.character(rep(1,item.count)))
}

cc.all.n.groups<-function(item.count=10)
{
  dissimilarity.matrix<-matrix(data = 1, nrow=item.count, ncol=item.count)
  diag(dissimilarity.matrix)<-NA
  pgroups<-PostHocGroups.dist(dissimiliarity.matrix = dissimilarity.matrix)
  best.row<-as.character(pgroups[order(-qual)][1,!colnames(pgroups) %in% c('id','sum.inner','sum.outer','n.inner','n.outer','qual'),with=FALSE])
  testthat::expect_equal(best.row, as.character(seq(1,item.count)))
}


test.simple.grouping<-function()
{
  dm<-matrix(c(
    NA, 0, 0, 0, 0, 0,
     0,NA, 0, 1, 1, 1,
     0, 0,NA, 1, 1, 1,
     0, 1, 1,NA, 0, 0,
     0, 1, 1, 0,NA, 0,
     0, 1, 1, 0, 0, NA),nrow=6)
  dm[dm==0]<- -1
  colnames(dm)<-c('1st','2nd','3rd','4th','5th','6th')
  pgroups<-PostHocGroups(as.dist(dm))
  ans<-get.table.struct(pgroups)
  car<-count.rows.and.columns(ans)
  tab<-get.excel.table(car)
  ans<-t(matrix(c( 1,2,
                   1,2,
                   1,3,
                   1,3,
                   1,3), nrow=2))
  testthat::expect_equal(bigmemory::as.matrix(tab$tab),ans)
}

test.recursive.grouping<-function()
{
  dm<-matrix(c(
  NA, 1, 0, 0, 0, 1, 1, 1, 1,
   1,NA, 1, 1, 1, 0, 0, 0, 0,
   1, 1,NA, 1, 1, 1, 1, 1, 1,
   1, 1, 1,NA, 1, 1, 1, 1, 1,
   1, 1, 1, 1,NA, 1, 1, 1, 1,
   1, 0, 1, 1, 1,NA, 1, 1, 1,
   1, 0, 1, 1, 1, 1,NA, 1, 1,
   1, 0, 1, 1, 1, 1, 1,NA, 1,
   1, 0, 1, 1, 1, 1, 1, 1,NA
  ),nrow=9)
  dm[dm==0]<- -1
  colnames(dm)<-LETTERS[1:9]
  pgroups<-PostHocGroups(as.dist(dm), max.recursive.level=2)
}


test<-function()
{
  nvar<-5
  classcnt<-2
  myprox<-matrix(sample(c(-2:2),nvar*nvar, replace = TRUE),nrow = nvar)
  colnames(myprox)<-1:nvar
  rownames(myprox)<-1:nvar
  solution<-sample(0:classcnt,nvar,replace=TRUE)
  class.to.ignore<-1
}
