## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(DMtest)
#load example data
data(beta)
dim(beta)
data("covariate")
dim(covariate)
#compute p-values 
out=dmvc(beta=beta,covariate=covariate)
head(out)

