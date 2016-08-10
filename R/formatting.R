#' Function that renders the the row into a more table-friendly struct, that can get relatively easly printed
#'
#' @return A list with two elements
#' \describe{
#'  \item{\code{common.items}}{character vector with values that will appear in rows
#'  on the left column in the result table
#'  }
#'  \item{\code{discrete.items}}{list with one item for each element on the right column(s).
#'  Each item is a list with two elements: \code{common.items} and \code{discrete.items} all described above
#'  }
#' }
#'
get.table.struct<-function(row)
{
  if (is.null(names(row)))
  {
    names<-colnames(row)
  } else {
    names<-names(row)
  }
  if ('_qual' %in% names)
  {
    names<-names[which(names!='_qual')]
    row<-row[,names]
  }
  #browser()
  zero.rows<-stringr::str_count(row,'^0$')
  discrete.rows<-stringr::str_count(row,'^[1-9].*$')
  ans<-list()
  if (sum(zero.rows)>0)
  {
    ans$common.items<-names[which(zero.rows==1)]
  }
  if (sum(discrete.rows)>0)
  {
    u.discrete.rows<-unique(stringr::str_match(row[which(discrete.rows==1)],'^(\\d+)')[,2])
    item<-u.discrete.rows[[2]]
    f<-function(item){
      idx<-which(stringr::str_count(row, paste0('^',item,'.*$'))==1)
      newrow<-stringr::str_match(row[idx],paste0('^',item,'\\.?(.*)$'))[,2]
      newrow[newrow=='']<-0
      names(newrow)<-names[idx]
      return(get.table.struct(newrow))
    }
    ans$discrete.items<-plyr::alply(u.discrete.rows,1,f)
    attributes(ans$discrete.items)<-NULL
  }
  return(ans)
}

#' This function gets table format rendered by \code{get.table.struct} and converts it into even more table-friendly
#' data representation.
#'
#' The return value is a list of max 4 elements:
#' * total row.count - if we keep on item on one row, this is the max of all columns in the resulting table
#' * common.column.count - either 1 or 0, depending on whether the common item column is present
#' * common.items - item numbers of all items that go to the common group (if one exists)
#' * discrete.column.count - total number of all columns that contain exclusive groups. More than one, if the
#'   search was recursive.
#' * discrete.items - list where each element describes every distinct exclusive group. The data format
#'   is the same as the main object (this is a recursively defined list)
count.rows.and.columns<-function(ans)
{
  #browser()
  comm.count<-0
  common.column.count<-0
  if (!is.null(ans$common.items))
  {
    comm.count<-length(ans$common.items)
    common.column.count<-1
  }
  discrete.items<-list()
  discr.row.count<-0
  discrete.column.count<-0
  subitem.has.common.column<-0
  for (i in seq_along(ans$discrete.items))
  {
    item<-ans$discrete.items[[i]]
    tmp<-count.rows.and.columns(item)
    if (tmp$common.column.count>0)
    {
      subitem.has.common.column<-1
    }
    discr.row.count<-discr.row.count+tmp$row.count
    discrete.column.count<-max(discrete.column.count,tmp$discrete.column.count)
    discrete.items[[i]]<-tmp
  }
  for (i in seq_along(discrete.items))
  {
    discrete.items[[i]]$common.column.count<-subitem.has.common.column
    discrete.items[[i]]$discrete.column.count<-discrete.column.count
  }
  return(list(
    row.count=max(comm.count, discr.row.count),
    common.column.count=common.column.count,
    discrete.column.count=discrete.column.count+subitem.has.common.column,
    common.items=ans$common.items,
    discrete.items=discrete.items
  ))
}

#' This is the final function that outputs the table in an Excel-friendly format.
#'
#' All cells that are marked with the same non-NA symbol, should be merged into one group.
get.excel.table<-function(car, legend=NULL, tab=NULL, coloffset=1, rowoffset=1, nextid=1)
{
#  browser()
  if (is.null(tab))
  {
    tab<-bigmemory::big.matrix(
      nrow = car$row.count,
      ncol = car$common.column.count+car$discrete.column.count,
      type='integer')
  }
  if (is.null(legend))
  {
    legend<-list() # id -> vector of classes
  }

  nextcol<-coloffset
  if (car$common.column.count>0)
  {
  #  browser()
    for (i in seq(from = rowoffset,length.out =  car$row.count))
    {
      tab[i,nextcol]<-as.integer(nextid)
    }
    coloffset<-coloffset+1
    legend[[as.character(nextid)]]<-car$common.items
    nextid<-nextid+1
  }
  cumsum<-0
  for (i in seq_along(car$discrete.items))
  {
    item<-car$discrete.items[[i]]
    ans<-get.excel.table(   car = item,
                            legend = legend,
                            tab = tab,
                            coloffset = coloffset,
                            rowoffset = rowoffset+cumsum,
                            nextid = nextid)
    legend<-ans$legend
    tab<-ans$tab
    nextid<-ans$nextid
    cumsum<-cumsum+item$row.count
  }
  return(list(tab=tab, nextid=nextid, legend=legend))
}

#' @export
print.PostHocGroup<-function(row)
{
  ans<-get.table.struct(row)
  car<-count.rows.and.columns(ans)
  tab<-get.excel.table(car)
  names<-character(0)
  if (car$common.column.count>0)
  {
    names<-'common items'
  }
  if (car$discrete.column.count>0)
  {
    names<-c(names, 'partitioned items')
  }
  ans<-bigmemory::as.matrix(tab$tab)
  if (is.null(dim(ans)))
  {
    dim(ans)<-c(length(ans),1)
  }
  colnames(ans)<-names
  rownames(ans)<-seq(1, nrow(ans))

  print(ans)
  cat("Legend:\n")
  for(l in seq_along(tab$legend))
  {
    cat(paste0("'",l,"' - ",paste0(tab$legend[[l]], collapse = ', ' ), '\n'))
  }


}
#' @export
latex.PostHocGroup<-function(obj,title=digest::digest(obj), file=NULL)
{
  ans<-get.table.struct(obj)
  car<-count.rows.and.columns(ans)
  tab<-get.excel.table(car)
  m<-ans<-bigmemory::as.matrix(tab$tab)
  if (is.null(dim(m)))
  {
    dim(m)<-c(length(m),1)
  }
  col.specs<-paste0(rep('l',ncol(m)),collapse = '|')
  latex_body<-rep('',nrow(m)+length(tab$legend)+3)
  latex_body_i<-3
  latex_body[[1]]<-paste0('\\newcommand{\\specialcellitem}[3][c]{%\n',
                          '  \\renewcommand{\\arraystretch}{#2}\\begin{tabular}[#1]{@{}c@{}}#3\\end{tabular}}\n')
  latex_body[[2]]<-paste0('\\begin{tabulary}{\\textwidth}{',col.specs,'}\n',
                          '    \\toprule\n',
                          paste0('    \\multicolumn{',ncol(m),'}{l}{Variable groups}\\\\\n'),
                          '    \\midrule\n')

  latex_headers<-c('multirow','booktabs', 'tabulary')
  calc.n.in.table<-function(x,item)
  {
    return(sum(m[seq(1,nrow(m)),x]==item))
  }
  new.group.pos<-function(y)
  {
    if (y==1)
      return(0)
    vec<-m[y,]!=m[y-1,]
    ans<-which.max(vec)
    if (!vec[[ans]])
      return(0)
    return(ans)
  }
  items.done<-rep(0,length(tab$legend))
  items.tex<-matrix('',nrow(m),ncol(m))
  matrix.elems.done<-matrix(FALSE,nrow(m),ncol(m))
  for(y in seq(1, nrow(m)))
  {
    for(x in seq(1, ncol(m)))
    {
      if (!matrix.elems.done[[y,x]])
      {
        my.group<-m[y,x]
        n.in.legend<-length(tab$legend[[my.group]])
        n.in.table<-calc.n.in.table(x,my.group)
        if (n.in.table != n.in.legend)
        {
          items.tex[[y,x]]<-paste0('\\multirow{',
                                    n.in.table,
                                    '}{*}{\\specialcellitem{',
                                    format(n.in.table/n.in.legend,digits=4, scientific=FALSE),
                                    '}{',
                                    paste0(tab$legend[[my.group]],collapse='\\\\'),
                                    '}}'
                      )
          matrix.elems.done[seq(y,y+n.in.table-1),x]<-TRUE
        } else {
          items.tex[[y,x]]<-tab$legend[[my.group]][[items.done[[my.group]]+1 ]]
          matrix.elems.done[[y,x]]<-TRUE
          items.done[[my.group]]<-items.done[[my.group]]+1
        }
      }

    }
  }
  for(y in seq(1, nrow(m)))
  {
    newcolpos<-new.group.pos(y)
    if (newcolpos>0)
    {
      latex_body[[latex_body_i]]<-paste0('    \\cline{',newcolpos,'-',ncol(m),'}\n')
      latex_body_i<-latex_body_i+1
    }
    latex_body[[latex_body_i]]<-paste0('    ',
                                       paste0(items.tex[y,], collapse = ' & '),
                                       ' \\\\\n    ')
    latex_body_i<-latex_body_i+1
  }

  latex_body[[latex_body_i]]<-paste0('    \\bottomrule\n',
                                     '\\end{tabulary}\n')
  latex_body_i<-latex_body_i+1
  if (is.null(file))
  {
    file=paste0(title, ".tex")
  }
  cat(paste0(latex_body,collapse = ''),file=file)
  return(structure(
    list(file=file, style=latex_headers),class='latex'))
}


# library(Hmisc)
# a<-readRDS('data/ex1.rds')
# obj<-GroupPostHocs(a, nice.labels = c("Date of birth", "Gender", "Ethnicity", "Local of Birth", "Father´s local of Birth", "Date of Diagnosis", "Date of First Symptoms"))
# f<-latex(obj)
# ll<-readLines(f$file)
# cat(paste0(ll,collapse = '\n'))
