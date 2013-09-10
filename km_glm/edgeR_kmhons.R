# This gets the analysis-specific options into the main glm multifactor script
modelmatrix <- model.matrix(
  ~ GrowthCondition + Treatment,# + GrowthCondition * Treatment,
  data=keyfile
  )

contrasts <- makeContrasts(
  "Suf-Exc"=GrowthConditionSufficient-Intercept,
  "Suf-Flu"=GrowthConditionSufficient-GrowthConditionFluctuating,
  "Flu-Exc"=GrowthConditionFluctuating-Intercept,
  levels=modelmatrix
  )
