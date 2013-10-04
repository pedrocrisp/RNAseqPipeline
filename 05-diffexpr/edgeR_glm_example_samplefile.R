# This gets the analysis-specific options into the main glm multifactor script

modelmatrix <- model.matrix(
  ~ 0 + Group,
  data=keyfile
  )

contrasts <- makeContrasts(
  "A-B"=A-B,
  levels=modelmatrix
)

FDR.cutoff <- 0.05
