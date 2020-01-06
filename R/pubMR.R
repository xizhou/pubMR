#' @slot name
#' @slot host
#' @slot db
#' @slot retmax
#' @slot rettype
#' @slot retmode
#' @slot email
#' @slot tool
#' @exportClass eutils
setClass("eutils", 
   representation=list(name="character",host="character",db="character",retmax="character", 
   rettype = "character",retmode="character",email="character",tool="character"), 
   prototype = list(db = "pubmed", retmode = "xml", tool = "reutils"))
#' @slot term
#' @exportClass search
setClass("search",contains="eutils",representation=list(term="character"), 
   prototype=list(host="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
   name="esearch",rettype ="count",retmax ="100"))
#' @slot id
#' @exportClass fetch
setClass("fetch",contains="eutils",representation=list(id="character"), 
   prototype=list(host = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"))
#' @export
setGeneric("Extractor", function(obj,...) standardGeneric("Extractor"))
#' @export
setMethod("Extractor", "search", function(obj,...)
{
   if (is.null(obj@term))
   {
      stop("No query term provided")
   }
   if (length(obj@term) > 1L)
   {
      obj@term <- paste(obj@term,collapse=" OR ")
   }
   field <- list(term=obj@term,db=obj@db,retmax=obj@retmax,rettype=obj@rettype, 
      retmode=obj@retmode,email=obj@email,tool=obj@tool)
   opts <- list(connecttimeout = 10)
   hg <- RCurl::basicHeaderGatherer()
   opts$headerfunction <- hg$update
   tg <- RCurl::basicTextGatherer()
   opts$writefunction <- tg$update
   if(nchar(obj@term) < 100)
   {
      url.field <- paste(RCurl::curlEscape(names(field)),RCurl::curlEscape(field),sep="=", 
                       collapse="&")
      url <- paste(obj@host, url.field, sep = "?")
      e <- RCurl::getURLContent(url = url, .opts = opts)
   }
   else e <- RCurl::postForm(obj@host,.params=field,.opts=opts,style="POST")
   res.txt <- as.character(tg$value())
   res.xml <- XML::xmlTreeParse(res.txt,asText=TRUE,useInternalNodes=T) 
   if(obj@rettype == "count")
       y <- XML::xpathSApply(res.xml,"/eSearchResult/Count", XML::xmlValue)
   else y <- XML::xpathSApply(res.xml, "//Id", XML::xmlValue)
   y
})
#' @export
setMethod("Extractor","fetch",function(obj,...)
{
   field <- list(id = paste(obj@id,collapse = ","),db=obj@db,retmax=obj@retmax, 
      rettype=obj@rettype,retmode=obj@retmode,email=obj@email,tool=obj@tool)
   opts <- list(connecttimeout=1000)
   hg <- RCurl::basicHeaderGatherer()
   opts$headerfunction <- hg$update
   tg <- RCurl::basicTextGatherer()
   opts$writefunction <- tg$update 
   if(length(obj@id) <= 200)
   {
      url.field <- paste(RCurl::curlEscape(names(field)),RCurl::curlEscape(field),
         sep="=",collapse="&")
      url <- paste(obj@host,url.field,sep="?")
      e <- RCurl::getURLContent(url=url,.opts=opts)
   }   
   else e <- RCurl::postForm(obj@host,.params=field,.opts=opts,style="POST")
   res.txt <- as.character(tg$value())
   res.xml <- XML::xmlTreeParse(res.txt,asText=TRUE,useInternalNodes=T,encoding="UTF-8") 
   res.xml
})

#' @export
AB <- function(query,input=NULL,output="ABprofile")
{
   format <- match.arg(output, c("ABprofile", "pmidlist", "xml"))
   if(is.null(input))
   {
      s <- new("search", term = query)
      #cat("we only extract English articles\n")
      s@retmax <- Extractor(s)
      s@rettype<- "uilist"
      pmidlist <- Extractor(s)
      if(format == "pmidlist"){
         return(pmidlist)
      }else{
      f <- new("fetch", id = pmidlist, retmax = as.character(length(pmidlist)))
      y <- Extractor(f)
   }}
   else if(length(grep(".xml",input))>0)
   {
      #inpout must be XML file
      y <- input
      y <- XML::xmlTreeParse(y,useInternalNodes = T, encoding = "UTF-8")
   }
   else
   {
      pmidlist <- input
      f <- new("fetch", id=pmidlist,retmax = as.character(length(pmidlist)))
      y <- Extractor(f)
      
   }
   if(format == "ABprofile")
        y <- parseAB(y)
   y
}

#' @slot PMID
#' @slot TI
#' @slot AB
#' @slot TA
#' @slot PDAT
#' @slot ISSN
#' @slot MH
#' @slot SH
#' @slot MAJR
#' @slot AU
#' @slot MS
#' @exportClass ABprofile
setClass("ABprofile", representation(PMID="character",TI="character",AB="character", 
   TA="character",PDAT="numeric",ISSN="character", 
   MH="list",SH="list",MAJR="list",AU="list",MS="list"))
#' @export
setMethod("show", "ABprofile",
function(object)
{
   cat("An object of class \"",class(object),"\"",
    " containing ",length(object@PMID),
    " articles", "\n",sep="")
   cat("with slot names: ",
      paste(slotNames(object),collapse=","),".","\n",sep="")
        
})

#' @export
setMethod("[", c("ABprofile", "numeric", "missing"),
    function(x, i, j, drop=TRUE)
{
   initialize(x,PMID=x@PMID[i],TI=x@TI[i],AB=x@AB[i],PDAT=x@PDAT[i],
      ISSN=x@ISSN[i],MH=x@MH[i],SH=x@SH[i],MAJR=x@MAJR[i],AU=x@AU[i],MS=x@MS[i])
})
#' @export
setMethod("[", c("ABprofile", "character"),
    function(x, i, drop=TRUE)
{
   i <- which(x@PMID %in% i)
   initialize(x,PMID=x@PMID[i],TI=x@TI[i],AB=x@AB[i],PDAT=x@PDAT[i],
      ISSN=x@ISSN[i],MH=x@MH[i],SH=x@SH[i],MAJR=x@MAJR[i],AU=x@AU[i],MS=x@MS[i])
})
#' @export
setMethod("$", "ABprofile",
    function(x, name)
{
    ## 'name' is a character(1)
    slot(x, name)
})

#' @export
parseAB <- function(doc)
{
   parse <- function(xmlnodeset, path)
   {    
      sapply(xmlnodeset, function(x)
      {
         temp = XML::xpathSApply(x, path, XML::xmlValue)
         if(length(temp)==0)
            y <- "NO Result Found"
         else y <- temp
         y
      })
   }  
   ArticleList <- XML::getNodeSet(doc, "//PubmedArticle")
   PmidList <- parse(ArticleList, ".//MedlineCitation/PMID[1]")
   TitleList <- parse(ArticleList, ".//ArticleTitle")
   AbstractList <- sapply(parse(ArticleList,".//AbstractText"), 
      function(x) paste0(x,collapse=" "))
   JournalList <- parse(ArticleList, ".//ISOAbbreviation")
   ISSN <- parse(ArticleList, ".//MedlineJournalInfo/ISSNLinking")
   YearList <- parse(ArticleList, ".//PubDate")
   YearList <- sapply(YearList, function(x)
   {
      if(x != "NO Result Found")
         x <- stringr::str_sub(x,1, 4)
      x
   })
   names(YearList) <- NULL
   MeshList <- parse(ArticleList, ".//MeshHeading//DescriptorName")
   SubmeshList <- sapply(parse(ArticleList, ".//QualifierName"), unique)
   MajrList <- sapply(ArticleList, function(x)
   {
       tmp <- XML::getNodeSet(x, 
          ".//DescriptorName[@MajorTopicYN='Y'] | .//QualifierName[@MajorTopicYN='Y']")
       if(length(tmp) == 0)
          y <- "NO Result Found"
       else y <- sapply(unique(sapply(tmp, XML::xmlParent)), 
          function(x) XML::xpathSApply(x, ".//DescriptorName", XML::xmlValue))
       y
   })
   AuthorList <- parse(ArticleList, 
      ".//AuthorList/Author/LastName | .//AuthorList/Author/ForeName")
   AuthorList <- sapply(AuthorList,function(x)
   {
      if(x[1] != "NO Result Found")
           x <- paste(x[(1:length(x)) %% 2 == 1], x[(1:length(x)) %% 2 == 0], sep = " ")
      x
   })
   MS <- sapply(ArticleList, function(x)
   {
       tmp1 <- XML::getNodeSet(x, ".//DescriptorName[@MajorTopicYN='Y']")
       tmp2 <- XML::getNodeSet(x, ".//QualifierName[@MajorTopicYN='Y']")
       y <- NULL
       if(length(tmp1)>0)
       {
          y1 <- sapply(tmp1, XML::xmlValue)
          y2 <- sapply(sapply(tmp1, XML::xmlParent), function(x) 
             XML::xpathSApply(x, ".//QualifierName", XML::xmlValue))
          res <- mapply(function(u,v)
             {
                if(length(v)) 
                   paste(u,v,sep="/")
                else u
             },u=y1,v=y2)
          res <- unlist(res)
          names(res) <- NULL
          y <- c(y,res)
      }
     if(length(tmp2)>0)
     {
        y1 <- sapply(sapply(tmp2, XML::xmlParent), 
           function(x) XML::xpathSApply(x, ".//DescriptorName", XML::xmlValue))
        y2 <- sapply(tmp2, XML::xmlValue)
        y <- c(y,paste(y1,y2,sep="/"))
     }
     y})                                   
    profile <- new("ABprofile", PMID = PmidList, TI = TitleList, AB = AbstractList, PDAT = as.numeric(YearList),
                   TA = JournalList, ISSN = ISSN, MH = MeshList, MAJR = MajrList, SH = SubmeshList, 
                   AU = AuthorList, MS=MS)

}

#' @export
pubtator <- function(pmidlist)
{
   .Extractor <- function(res.txt)
   {
       x <- try(XML::xmlTreeParse(res.txt,asText=TRUE,useInternalNodes=T),silent=TRUE)
       if(class(x)[1L]=="try-error")
          y <- NA
       else if(length(XML:::getNodeSet(x, "//annotation")) == 0)
          y <- NA
       else
       {
          tmp <- XML::getNodeSet(x,"//id")
          pmid <- XML::xmlSApply(tmp,XML::xmlValue)
          ann <- XML::getNodeSet(x,"//annotation")
          id <- name <- type <- rep(NA,length(ann))
          for(i in seq(ann))
          {
             tmp1 <- XML::getNodeSet(ann[[i]], ".//infon[@key='identifier']|.//infon[@key='Identifier']")
             id[i] <- ifelse(length(tmp1)==0,NA,XML::xmlValue(tmp1[[1]]))
             tmp2 <- XML::getNodeSet(ann[[i]], ".//text")
             name[i] <- ifelse(length(tmp2)==0,NA,XML::xmlValue(tmp2[[1]]))
             tmp3 <- XML::getNodeSet(ann[[i]], ".//infon[@key='type']")
             type[i] <- ifelse(length(tmp3)==0,NA,XML::xmlValue(tmp3[[1]]))
          }
             y <- data.frame(pmid,type,name,id,stringsAsFactors=F)
       }  
       y   
   }
   PubtatorAPI <- function(pmidlist)
   {
      opts <- list(connecttimeout=1000)
      hg <- RCurl::basicHeaderGatherer()
      opts$headerfunction <- hg$update
      tg <- RCurl::basicTextGatherer()
      opts$writefunction <- tg$update
      e <- RCurl::getURLContent(url=
          paste0("https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocxml?pmids=", 
          paste(pmidlist, collapse = ",")), .opts = opts)
      res.txt <- as.character(tg$value())
      res.txt <- unlist(str_split(res.txt,"<document>"))
      res.txt <- paste0("<document>",res.txt[-1])
      res.txt[length(res.txt)] <- stringr::str_replace(res.txt[length(res.txt)],"</collection>", "")   
      annota.list <- lapply(res.txt,.Extractor)
      Sys.sleep(0.01)
      return(annota.list)
   } 
   if(length(pmidlist) <= 100)
       res <- PubtatorAPI(pmidlist)
   else
   {
      res <- list()
      s <- split(pmidlist, ceiling(seq_along(pmidlist)/100))
      res <- lapply(s,PubtatorAPI)
      res <- unlist(res,recursive=FALSE)
   }
   res <- res[!is.na(res)]
   res <- data.table::rbindlist(res)
   res
}