## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## ----install------------------------------------------------------------------
#  devtools::install_github("FranjoIM/PackageName", ref = "v0.1.0", build_vignettes = TRUE)

## ----setup--------------------------------------------------------------------
#  library(GWASr)

## ----vignette-----------------------------------------------------------------
#  vignette("GWASr", package = "GWASr")

## ----Genomic lambda-----------------------------------------------------------
#  GLambda(ADHDMetaP$P)

## ----Genomic lambda tes-------------------------------------------------------
#  GLambda(0.5)
#  #> 1

## ----Geenomic lambda arguments------------------------------------------------
#  ?GLambda

## ----Manhattan plot-----------------------------------------------------------
#  plot_manhattan(ADHDMeta)

## ----Manhattan plot arguments-------------------------------------------------
#  ?plot_manhattan

## ----QQ plot------------------------------------------------------------------
#  plot_qq(DF)

## ----QQ plot arguments--------------------------------------------------------
#  ?plot_qq

## ----ADHDMeta-----------------------------------------------------------------
#  ADHDMeta
#  names(ADHDMeta)
#  #> [1] "CHR" "SNP" "BP"  "OR"  "P"

## ----ADHDMetaP----------------------------------------------------------------
#  ADHDMetaP
#  names(ADHDMetaP)
#  #> [1] "SNP" "P"

