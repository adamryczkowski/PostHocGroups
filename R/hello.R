# This package provides means for simple grouping of lists of things,
# for which we can define a relation 'is.different(thing1, thing2)'
#
# The package internally creates annealed dissimilarity matrix, which is then
# basis for clustering.
#
# Author: Adam Ryczkowski
###############################################################################

if (!require("R.utils")) install.packages("R.utils")
R.utils::use("plyr")

#public:
group.results<-function(co,fkcja01)
{
  dobreidx<-which(laply(co,class)=='fitdist') #indeksy niezerowych
  dobreco<-co[dobreidx]
  m<-symmetric.list.outer(dobreco,fkcja01) #macierz symetryczna podobieństwa obiektów
  res<-znajdzrozwiazanie(m) #ale jeszcze trzeba to odnieść do niezerowych
  out<-rep(NA,length(co))
  out[dobreidx]<-res$part
  list(part=out,qual=res$qual)
}

makeprox<-function(mycutree){
  #Funkcja tworzy macierz odległości opartą na wyniku cutree klastrowania
  #wynik?w post-hoc
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

util.f.inner.group<-function(x)
{
  sqrt(abs(x))
}


qualofpartition<-function(mycutree, myprox,
                          prog.alfa=0.05,
                          util.f.inner.group = base::identity,
                          util.f.extra.group=base::identity)
{
  #Funkcja zwraca jako?? danego pogrupowania.
  #mycutree<-mycutrees[38,]
#  browser()
  otherprox<-as.matrix(makeprox(as.integer(mycutree)))
#  myprox.m<-as.matrix(myprox)
  myprox.m<-myprox
  filtr<-otherprox==0 & lower.tri(myprox.m)
  m0<-util.f.inner.group(myprox.m[filtr]-log(prog.alfa))
  m0<-mean(m0)

  filtr<-otherprox==1 & lower.tri(myprox.m)
  m1<-util.f.inner.group(myprox.m[filtr]-log(prog.alfa))
  m1<-mean(m1)

  ans<-mean(c(m0,m1), na.rm=TRUE)
  return(ans)
}


znajdzrozwiazanie<-function(myprox,means)
{
  myprox<-as.dist(m=myprox)
  bestqual<-1 #najgorsza możliwa
  sig.count<-plyr::aaply(as.matrix(myprox),2,sum) #wektor z liczbą istotnotnych różnic dla każdego z elementów
  mylen<-attr(myprox,"Size")
  ord<-order(means)
  smile<-sig.count[ord] #Items in the middle should have half the
  #amount of significant differences, than variables on the edges, that's why the name.
  #Now we will correct for the effect:

  f<-function(x)
  { #Function with "ideal" smile: f(0)=1, f(mylen-1)=1, f((mylen-1)/2)=0.5
    return((2*x*x)/((mylen-1)*(mylen-1))-(2*x)/(mylen-1)+1)
  }
  x<-seq(0,mylen-1)
  A<-coef(nls(smile ~ A*f(x), start=list(A=mean(smile)))) #fitted the only degree of freedom: amplitude of the function.
  smile2<-plyr::aaply(x,1,function(x){f(x)*A}) #Idealistic, fitted smile
  smile.cor<-smile-smile2 #Corrected smile.
  sig.count.cor<-smile.cor[Matrix::invPerm(ord)]

  #Empty data.table to store results
  solutions<-as.data.table(read.table(text="",colClasses = 'integer',col.names = 1:mylen))
  colnames(solutions)<-as.character(c(1:mylen))
  for(i in seq(2, floor(sqrt(mylen)))) #Liczba klastrów, na które dzielimy wszystkie zmienne.
    #Liczba klastrów determinuje proporcję podziału zmiennych na 2 grupy.
  {
    n<-kmeans(sig.count, i)
    ocenters<-order(n$centers)
    for (j in seq(1, length(n$centers)-1)) # j = pozycja odcięcia między klastrami zmiennych na te,
      # które zostają a zmiennymi,
      # które partycjonujemy. j = pozycja początku zakresu do partycjonowania; j=1 -> brak zmiennych, które zostają.
    {
 #     cat(paste('j=',j))
      whichc.nr<-ocenters[seq(j,length(ocenters))] #Numery grup grup, które chcemy zostawić dla algorytmu
      whichc<-as.numeric(names(n$cluster))[n$cluster %in% whichc.nr] #Które elementy (columny) zostawiamy dla algorytmu?
      if (length(whichc)>1)
      {
        myprox2<-as.dist(as.matrix(myprox)[whichc,whichc])
        mylen2<-attr(myprox2,"Size")
        wyn<-hclust(d=myprox2)
        solutions<-findbestpartition(wyn,myprox2,solutions)
      }
    }
#    cat(paste0('i=',i))
  }
#  browser()
  solutions<-unique(solutions)

  #TODO: Add column that calculates quality
}

#Funkcja zwraca 0, gdy dwa zakresy si? zaz?biaj? i 1, gdy zakresy nie maj? punkt?w wsp?lnych
interval.overlaps<-function(para1,para2)
{
  if(sum(is.na(c(para1,para2)))>0)
    return(NA)

  if(para1[1]>para2[2])
    return(1)

  if(para2[1]>para1[2])
    return(1)

  if(para1[2]<para2[1])
    return(1)

  if(para2[2]<para1[1])
    return(1)

  return(0)
}



