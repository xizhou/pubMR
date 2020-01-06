pubMR
==========
pubMR是R平台下一个高效的PubMed文本挖掘工具，集合了：检索下载、解析抽取、基本统计、多维矩阵、论文相似、热点分析、概念识别和网络分析等多种功能。

pubMR is an R package designed for text mining of the PubMed database. Additionally, it provide some highly customized metics to evaluate and visualize results for downstream analysis.

## Installation

```r
devtools::install_github("xizhou/pubMR")
```

## Vignettes
[pubMR.pdf](https://github.com/xizhou/pubMR/tree/master/vignettes/pubMR.pdf)

## Usage
A quick start:
```r
library(pubMR)
m <- '"neoplasms"[MeSH Terms] AND "serine/metabolism"[Mesh Terms] AND ("2017/01/01"[PDAT] : "2018/12/31"[PDAT])'
obj <- AB(query=m,output='ABprofile')
p <- pubtator(obj@PMID)
```
Import an xml file downloaded from the PubMed database into R program:

<p align="center">
  <img src="https://github.com/xizhou/pubMR/blob/master/screenshot.png?raw=true" alt="xml"/>
</p>

```r
library(pubMR)
obj <- AB(input="pubmed_result.xml")
```