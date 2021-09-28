## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----download-----------------------------------------------------------------
#CellMarker
cellmarker = read.table(url("http://biocc.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt"), header=TRUE, sep="\t")

#CEllmatch
sysdata_cellmatch=tempfile()
curl::curl_download(url="https://github.com/ZJUFanLab/scCATCH/raw/master/R/sysdata.rda", sysdata_cellmatch)
load(sysdata_cellmatch)
ls()
cellmatch = CellMatch
unlink(sysdata_cellmatch)

## ----filter cellmarker--------------------------------------------------------

### Filtering the cellmarker object

## Remove NA values
cellmarker = cellmarker[-which(is.na(cellmarker$geneSymbol)), ]
head(cellmarker)

# Retrieve all the genes per cell population: each line corresponds to a particular cell populations. The gene markers are separated by a coma.
allgenes = lapply(sapply(cellmarker$geneSymbol,function(x){strsplit(as.character(x),", ")}), function(x){genes = unique(gsub("\\[|]", "", x)) ; return(genes)})

# Prepare the names of cell population mentionning the tissue type and cell type
cellID = paste(gsub(" ", "_", cellmarker$tissueType),  gsub(" ", "_", cellmarker$cellType), gsub("[+]", "plus", gsub("[(]|[)]|-","", gsub(" ", "_", cellmarker$cellName))), sep=".")

# Transform the df into a list. Note that some cell populations are specified several times in the file, with common genes. In addition we remove the NA values.
liter_cellmarker = tapply(cellmarker$geneSymbol, cellID, function(x){genes = unique(unlist(strsplit(as.character(x),", "))) ; if(sum(genes=="NA")>0){genes = genes[-which(genes=="NA")]} ; return(genes)})

## Filtering for the number of genes.

nbgenes = unlist(lapply(liter_cellmarker, length))
liter_cellmarker_filter = liter_cellmarker[which(nbgenes>=3 & nbgenes<=200)]

length(liter_cellmarker_filter)
liter_cellmarker_filter[1:3]


## ----filter cellmatch---------------------------------------------------------

### Filtering the cellmarker object

## Keep only human markers
cellmatch = cellmatch[cellmatch$speciesType == "Human",]
head(cellmatch)

## Prepare the names of cell population mentionning the tissue type and cell type
cellID = paste(gsub(" ", "_", cellmatch$tissueType),  gsub(" ", "_", cellmatch$cellType), gsub("[+]", "plus", gsub("[(]|[)]|-","", gsub(" ", "_", cellmatch$cellName))), sep=".")

## Transform the df into a list. 
liter_cellmatch = tapply(cellmatch$geneSymbol, cellID, list)
liter_cellmatch = lapply(liter_cellmatch, function(x) unique(x))

## Filtering for the number of genes.
nbgenes = unlist(lapply(liter_cellmatch, length))
liter_cellmatch_filter = liter_cellmatch[which(nbgenes>=3 & nbgenes<=200)]

length(liter_cellmatch_filter)
liter_cellmatch_filter[1:3]


## ----save, eval=FALSE---------------------------------------------------------
#  saveRDS(liter_cellmarker_filter, "cellmarkerDB.RDS")
#  
#  saveRDS(liter_cellmatch_filter, "cellmatchDB.RDS")
#  

