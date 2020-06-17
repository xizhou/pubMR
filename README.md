

pubMR
==========
pubMR是R平台下一个高效的PubMed文本挖掘工具，集合了：检索下载、解析抽取、基本统计、多维矩阵、论文相似、热点分析、概念识别和网络分析等多种功能。

pubMR is an R package designed for text mining of the PubMed database. Additionally, it provide some highly customized metics to evaluate and visualize results for downstream analysis.


## Authors

[崔雷] (Cuilei)

[周晓北] (Zhou Xiaobei)

## Installation
### for Linux and Mac:

```r
install.packages("devtools")
devtools::install_github("xizhou/pubMR")
```

### for WIN user:

```r
install.packages("devtools")
install.packages("RCurl", type="win.binary")
install.packages("stringi",type="win.binary")
devtools::install_github("xizhou/pubMR",upgrade="never",force=TRUE)
```

<!-- 
## Vignettes
pubMR.pdf](https://github.com/xizhou/pubMR/tree/master/vignettes/pubMR.pdf)
 -->
 
## Usage
### A quick start:
```r
library(pubMR)
m <- '"neoplasms"[MeSH Terms] AND "serine/metabolism"[Mesh Terms] AND ("2017/01/01"[PDAT] : "2018/12/31"[PDAT])'
obj <- txtList(input=m)
```

```text
obj
An object of class "txtList" containing 95 articles
with slot names: PMID,TI,AB,JT,DP,ISSN,MH,SH,MAJR,AU.
```
All text information are storing in format of **txtList**. **txtList** includes PMID, TI(title), AB(abstract), JT(journal title), DP(date publish), ISSN,
MH(MeSH heading), SH(MeSH subheading), MAJR(major topic, asterisks on MeSH headings and subheadings) and AU(author).

Any information (e.g., PMID) an be extracted like:
```text
obj@PMID
 [1] "30579699" "30527807" "30487607" "30413706" "30409761" "30217568"
 [7] "30209241" "30166347" "30157431" "30124210" "30118680" "30106121"
[13] "30051594" "30035852" "30028969" "30022161" "30017192" "29959989"
[19] "29949874" "29873416" "29702197" "29684847" "29672864" "29663362"
[25] "29636461" "29626472" "29615789" "29514256" "29507620" "29507618"
[31] "29452640" "29415992" "29408512" "29393406" "29361117" "29328365"
[37] "29301793" "29300860" "29225033" "29196224" "29180469" "29175460"
[43] "29133416" "29117943" "29038521" "29038347" "28960612" "28951458"
[49] "28945225" "28931725" "28919412" "28864417" "28787469" "28780319"
[55] "28765091" "28745319" "28708069" "28614715" "28607489" "28601431"
[61] "28555617" "28504189" "28500236" "28425994" "28422952" "28415597"
[67] "28412963" "28355568" "28350087" "28349780" "28328320" "28294115"
[73] "28262924" "28259896" "28236852" "28232235" "28192480" "28176663"
[79] "28154322" "28137613" "28130399" "28090628" "28077578" "28038466"
[85] "27996159" "27959531" "27939839" "27927748" "27902483" "27890797"
[91] "27890529" "27889207" "27777073" "27297361" "27748765"
```
 
### Advanced operations:
**Load/save files:** 
- Import a "PubMed" file downloaded from the PubMed database into R program:

<p align="center">
  <img src="https://github.com/xizhou/pubMR/blob/master/screenshot.png?raw=true" alt="PubMed"/>
</p>

```r
library(pubMR)
obj <- txtList(input="pubmed-neoplasmsM-set.txt",inputType="PubMed")
```
- import an "xml" file:
```r
library(pubMR)
obj <- txtList(input="pubmed_result.xml",inputType="xml")
```
- Save as an "xml" file:
```r
library(pubMR)
library(XML)
m <- '"neoplasms"[MeSH Terms] AND "serine/metabolism"[Mesh Terms] AND ("2017/01/01"[PDAT] : "2018/12/31"[PDAT])'
obj1 <- txtList(input=m,outputType='xml')
saveXML(obj1,file="result.xml")
```
**Produce co-occurrence-matrices:** 
- A MAJR-PMID co-occurrence-matrix can be generated by:
```r
library(data.table)
library(tidyr)
library(pubMR)
m <- '"neoplasms"[MeSH Terms] AND "serine/metabolism"[Mesh Terms] AND ("2017/01/01"[PDAT] : "2018/12/31"[PDAT])'
obj <- txtList(input=m)
obj1=data.table(PMID=obj@PMID,MAJR=obj@MAJR)
MAJR <- obj1[,MAJR]
idx <- sapply(MAJR,is.null)
obj1 <- obj1[!idx,]
obj1 = obj1 %>% unnest(MAJR) 
V <- table(obj1[,c("MAJR","PMID")])
```
```text
V
                                                 PMID
MAJR                                                27297361 27748765 27777073 27889207 27890529
  A Kinase Anchor Proteins/metabolism                    0        0        0        0        0
  Acquired Immunodeficiency Syndrome/metabolism          0        0        0        0        0
  Activating Transcription Factor 4/metabolism           0        0        0        0        0
  Adaptor Proteins, Signal Transducing/genetics          0        0        0        0        0
  Adaptor Proteins, Signal Transducing/metabolism        0        0        0        0        0
```

- A MAJR-MAJR co-occurrence-matrix can be generated by:
```r
V1 <- crossprod(t(V))
```
```text
V1                                                 MAJR
MAJR                                                A Kinase Anchor Proteins/metabolism Acquired Immunodeficiency Syndrome/metabolism Activating Transcription Factor 4/metabolism
  A Kinase Anchor Proteins/metabolism                                               1                                             0                                            0
  Acquired Immunodeficiency Syndrome/metabolism                                     0                                             1                                            0
  Activating Transcription Factor 4/metabolism                                      0                                             0                                            1
  Adaptor Proteins, Signal Transducing/genetics                                     0                                             0                                            0
  Adaptor Proteins, Signal Transducing/metabolism                                   0                                             0                                            0
```
- A Disease-Chemical-MAJR co-occurrence-matrix can be generated by:
```R
meshtree <- "https://github.com/xizhou/pubMR/raw/master/meshtree2019.Rdata"
load(url(meshtree))
```
```R
m <- '"neoplasms"[MeSH Terms] AND "serine/metabolism"[Mesh Terms] AND ("2017/01/01"[PDAT] : "2018/12/31"[PDAT])'
obj <- txtList(input=m)
obj1=data.table(PMID=obj@PMID,MAJR=obj@MAJR)
MAJR <- obj1[,MAJR]
idx <- sapply(MAJR,is.null)
obj1 <- obj1[!idx,]
obj1 = obj1 %>% unnest(MAJR) 
V <- table(obj1[,c("MAJR","PMID")])
V1 <- crossprod(t(V))
```
```R
nms <- rownames(V1)
nms <- gsub("\\/.*","",nms)
idr <- which(nms %in% meshtree[class=="D",mesh])
idc <- which(nms %in% meshtree[class=="C",mesh])
V2 <- V1[idr,idc]
```
**Combine with other tools:**
- The co-occurrence-matrix can be put into "gCluto" program:
```r
library(slam)
x <- V1
x <- x[!rowSums(x)==0,] 
x <- x[!colSums(x)==0,] 
slam::write_stm_CLUTO(x,file="dat.mat")
```
<p align="center">
  <img src="https://github.com/xizhou/pubMR/blob/master/fig.png?raw=true" alt="gcluto"/>
</p>
