#Estimate and save the parameters of the nb
estAndSave = function(physeq,  phyVar, tag, map = "HPCbootstrap",...){
  paramsName = paste0("params", tag)
  assign(paramsName, estNBparamsPhylo(physeq, phyVar,...))
  save(list = paramsName, file = file.path(map, paste0("params/",tag, "params.RData")))
  return(invisible())
}
