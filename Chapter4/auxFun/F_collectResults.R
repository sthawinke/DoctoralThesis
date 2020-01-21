loadBootStats = function(tag, map = file.path("HPCbootstrap/bootResults", tag)){
  List = list.files(map)
  tmp = lapply(List, function(x){
    load(file.path(map,x));bootTestStats
  })
  names(tmp) = gsub(".RData","",gsub("boot","",List))
  tmp
}
loadObsStats = function(tag, filePath = file.path("HPCbootstrap/realResults", paste0(tag, ".RData"))){
  load(filePath)
  return(get(paste0("smoothTest", tag)))
}
loadTestStats = function(tag, map = file.path("HPCbootstrap/bootResults", tag), nCores = nCores){
  List = list.files(map)
  tmp = lapply(List, function(x){
    load(file.path(map,x));bootTestStats
  })
  names(tmp) = gsub(".RData","",gsub("boot","",List))
  tmp
}