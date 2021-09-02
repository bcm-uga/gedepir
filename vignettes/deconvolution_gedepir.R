## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----datas--------------------------------------------------------------------
D = readRDS(url("https://figshare.com/ndownloader/files/30587814"))
dim(D)
head(D)[1:5, 1:5]


## ----run_norm-----------------------------------------------------------------
library(gedepir)
 D_norm=run_norm(mix_matrix = D,method = "RPM")

## ----run_trans----------------------------------------------------------------
 D_trans= run_trans(mix_matrix = D_norm,method = "log2")

## ----run_featsel--------------------------------------------------------------
 D_fsel= run_featsel(mix_matrix = D_trans,method = "cv1000")

## ----deconv-------------------------------------------------------------------
results_NMF = run_deconv(mix_matrix = D_fsel,k= 9, method = "NMF")

## ----enrich-------------------------------------------------------------------
library(fgsea)
database= readRDS(system.file("data", "Tissue.CellTypes.markers.RDS", package = "gedepir"))
gedepir::enrich(results_NMF$T_matrix,pathways = database,ICAbased = FALSE, fdr = FALSE)

## ----enrich-------------------------------------------------------------------
library(fgsea)
database= readRDS(system.file("data", "Tissue.CellTypes.markers.RDS", package = "gedepir"))
gedepir::enrich(results_NMF$T_matrix,pathways = database,ICAbased = FALSE, fdr = FALSE)

## ----pipeline-----------------------------------------------------------------
library(magrittr)
results_ICA= D %>%
  run_norm(method = "RPM") %>%
  run_trans(method = "linear") %>%
  run_deconv(method = "ICA")
gedepir::enrich(results_ICA$T_matrix,database,ICAbased = TRUE)

