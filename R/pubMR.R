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
txtList <- function(input,inputType="query",outputType="txtList")
{
   outputType <- match.arg(outputType, c("txtList", "pmidList", "xml"))
   inputType <- match.arg(inputType, c("query", "pmidList", "xml","PubMed"))
   if(inputType=="query")
   {
      s <- new("search", term=input)
      s@retmax <- Extractor(s)
      s@rettype<- "uilist"
      pmidList <- Extractor(s)
      if(outputType=="pmidList")
      {
         y <- pmidList
      }
      else 
      {
         f <- new("fetch",id=pmidList,retmax=as.character(length(pmidList)))
         y <- Extractor(f)
         if(outputType=="txtList")
            y <- .txtList(y)
      }
   }
   else if(inputType=="PubMed")
   {
      y <- readPB(input)
      y <- cLToS4(y,"txtList")
   }
   else if(inputType=="pmidList")
   {
      pmidList <- input
      f <- new("fetch", id=pmidList,retmax=as.character(length(pmidList)))
      y <- Extractor(f) 
      if(outputType=="txtList")
         y <- .txtList(y)     
   } 
   
   else if(inputType=="xml")
   {
      #inpout must be XML file
      y <- XML::xmlTreeParse(input,useInternalNodes=T,encoding="UTF-8")
      if(outputType=="txtList")
         y <- .txtList(y)
   }
   y
}

#' @slot PMID
#' @slot TI
#' @slot AB
#' @slot JT
#' @slot DP
#' @slot ISSN
#' @slot MH
#' @slot SH
#' @slot MAJR
#' @slot AU
#' @exportClass txt
setClass("txtList", representation(PMID="character",TI="character",AB="character", 
   JT="character",DP="character",ISSN="character", 
   MH="list",SH="list",MAJR="list",AU="list"))
#' @export
setMethod("show", "txtList",
function(object)
{
   cat("An object of class \"",class(object),"\"",
    " containing ",length(object@PMID),
    " articles", "\n",sep="")
   cat("with slot names: ",
      paste(slotNames(object),collapse=","),".","\n",sep="")
        
})

#' @export
setMethod("[", c("txtList", "numeric", "missing"),
    function(x, i, j, drop=TRUE)
{
   initialize(x,PMID=x@PMID[i],TI=x@TI[i],AB=x@AB[i],DP=x@DP[i],JT=x@JT[i],
      ISSN=x@ISSN[i],MH=x@MH[i],SH=x@SH[i],MAJR=x@MAJR[i],AU=x@AU[i])
})
#' @export
setMethod("[", c("txtList", "character"),
    function(x, i, drop=TRUE)
{
   i <- which(x@PMID %in% i)
   initialize(x,PMID=x@PMID[i],TI=x@TI[i],AB=x@AB[i],DP=x@DP[i],JT=x@JT[i],
      ISSN=x@ISSN[i],MH=x@MH[i],SH=x@SH[i],MAJR=x@MAJR[i],AU=x@AU[i])
})
#' @export
setMethod("$", "txtList",
    function(x, name)
{
    ## 'name' is a character(1)
    slot(x, name)
})

#' @export
.txtList <- function(doc)
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
   JournalList <- parse(ArticleList, ".//Title")
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
   #MajrList <- sapply(ArticleList, function(x)
   #{
   #    tmp <- XML::getNodeSet(x, 
   #       ".//DescriptorName[@MajorTopicYN='Y'] | .//QualifierName[@MajorTopicYN='Y']")
   #    if(length(tmp) == 0)
   #       y <- "NO Result Found"
   #    else y <- sapply(unique(sapply(tmp, XML::xmlParent)), 
   #       function(x) XML::xpathSApply(x, ".//DescriptorName", XML::xmlValue))
   #    y
   #})
   AuthorList <- parse(ArticleList, 
      ".//AuthorList/Author/LastName | .//AuthorList/Author/ForeName")
   AuthorList <- sapply(AuthorList,function(x)
   {
      if(x[1] != "NO Result Found")
           x <- paste(x[(1:length(x)) %% 2 == 1], x[(1:length(x)) %% 2 == 0], sep = " ")
      x
   })
   MAJR <- sapply(ArticleList, function(x)
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
    new("txtList",PMID=PmidList,TI=TitleList,AB=AbstractList,
                   DP=YearList,JT=JournalList,ISSN=ISSN,MH=MeshList,
                   MAJR=MAJR,SH=SubmeshList,AU=AuthorList)
}

readPB <- function(x)
{ 
   require(data.table)
   f <- readLines(x)
   f <- as.data.table(f)
   names(f) <- "V1"
   f <- split(f,cumsum(f==""))
   f <- lapply(f,function(x) x[,paste(V1,collapse=""), cumsum(!grepl("      ",V1))])
   f <- lapply(f,function(x) x[,gsub("       "," ",V1)])
   y <- list()
   for(i in seq(f))
   {
      x <- f[[i]]
      id <- grep("FAU - ",x)
      AU <- gsub("FAU - ","",x[id])
      id <- grep("JT  - ",x)
      JT <- gsub("JT  - ","",x[id])
      id <- grep("PMID- ",x)
      PMID <- gsub("PMID- ","",x[id])
      id <- grep("MH  - ",x)
      mesh <- gsub("MH  - ","",x[id])
      SH <- strsplit(mesh,"/")
      SH <- lapply(SH,function(x) x[-1L])
      SH <- SH[sapply(SH,length)>0]
      SH <- unique(gsub("\\*","",unlist(SH)))
      MH <- strsplit(mesh,"/")
      MH <- lapply(MH,function(x) x[1L])
      MH <- unique(unlist(MH))
      MAJR <- mesh[grep("\\*",mesh)]
      MAJR <- strsplit(MAJR,"/")
      MAJR <- lapply(MAJR,function(x) 
                        {id <- grep("\\*",x)
                        if(any(id%in%1))
                        {
                            if(length(id)==1)
                               y <- x[1]
                            else
                               y <- paste(x[1],x[-1L],sep="/")
                        }
                        else y <- paste(x[1],x[id],sep="/")
                        y})
      MAJR <- unlist(MAJR)
      MAJR <- gsub("\\*","",MAJR)
      id <- grep("AB  - ",x)
      AB <- gsub("AB  - ","",x[id])
      id <- grep("TI  - ",x)
      TI <- gsub("TI  - ","",x[id])
      id <- grep("DP  - ",x)
      DP <- gsub("DP  - ","",x[id])
      DP <- gsub(" .*","",DP)
      id <- grepl("IS  - ",x,fixed=TRUE)&grepl(" (Linking)",x,fixed=TRUE)
      IS <- gsub("IS  - ","",x[id])
      ISSN <- gsub(" (Linking)","",IS,fixed=TRUE)
      tmp <- list(AU=AU,JT=JT,PMID=PMID,SH=SH,MH=MH,
         MAJR=MAJR,TI=TI,AB=AB,DP=DP,ISSN=ISSN)
      tmp[sapply(tmp,length)==0] <- NA
      y[[i]] <- tmp
   }
   z <- c("AU","JT","PMID","SH","MH","MAJR","TI","AB","DP","ISSN")
   names(z) <- z
   lapply(z,function(x) sapply(y,"[[",x))
}


cLToS4 <- function (list, class) 
{
##originally import from "spectralAnalysis" pkg
   S4Object <- new(class)
   slotNames <- slotNames(S4Object)
   listNames <- names(list)
   checkNames <- all(listNames %in% slotNames)
   checkEqualNameLength <- length(listNames) == length(slotNames)
   if (!checkNames)
   {
      stop("names of 'list' do not match 'class'")
   }
   if (!checkEqualNameLength)
   {
      warning("some elements will be lost by conversion to S4 object")
   }
   for (slot in slotNames)
   {
      slot(S4Object, slot, check = FALSE) <- list[[slot]]
   }
   isValidObject <- validObject(S4Object)
   return(S4Object)
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
      res.txt <- unlist(stringr::str_split(res.txt,"<document>"))
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