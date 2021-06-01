#' Function to test enrichment of biological terms on each estimated component (matrix T).


# function to reorient the components if ICA-based deconvolution algorithm was performed
.orient_funct <- function(S) {
  orient <-
    apply(S, 2, function(x) {
      if (min(x) < -3 & max(x) > 3) {
        ifelse (sum(x > 3)  < sum(x < -3), -1, 1)
      } else {
        ifelse (sum(x > 2)  < sum(x < -2), -1, 1)
      }
    })
  S <- as.matrix(S)  %*% diag(orient)
  return(S)
}


# function to run the fgsea function and return a dataframe of the results
.enrichfun = function(pathways, genescores, showCategory = 10, showLeadingGenes = FALSE, fdr = TRUE, multilevel = FALSE){
    genes=sort(genescores,decreasing=TRUE)
    if(multilevel == FALSE){
        fgseaRes = data.frame(fgsea(pathways, genes, minSize=2, maxSize=200))
    } else if(multilevel == TRUE){
        fgseaRes = data.frame(fgseaMultilevel(pathways, genes, minSize=2, maxSize=200, eps=0))
    }
    pv = ifelse(fdr == TRUE,"padj","pval")
    res = fgseaRes[,c("pathway", "padj","pval", "ES")]
    if(showLeadingGenes == TRUE){
        res$leadingEdge = unlist(lapply(fgseaRes$leadingEdge, paste, collapse = "_"))
    }
    res = res[which(res[,pv]<0.05),]
    res = res[order(res[,pv], decreasing = FALSE),]
    respos = res[res[,"ES"]>0,]
    if(nrow(respos)==0){
        p = print("No significant cell-type enrichment")
    } else {
        respos = respos[1:min(showCategory, nrow(respos)),]
        resf = respos
    }
}

enrich = function(mydata, pathways = NULL, showCategory = 10, showLeadingGenes = FALSE, fdr = TRUE, multilevel = FALSE, ICAbased = FALSE){
    if(ICAbased == TRUE){
        mydata = .orient_funct(mydata)
    }
    apply(mydata, 2, function(x){.enrichfun(pathways = pathways, genescores = x, showCategory = showCategory, showLeadingGenes = showLeadingGenes, fdr = fdr, multilevel = multilevel)})
}

