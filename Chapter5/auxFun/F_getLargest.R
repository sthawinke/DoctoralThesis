#Get the largest objects
getLargest = function(num = 30){
    sort(decreasing = TRUE, sapply(ls(envir = globalenv()), function(x) object.size(get(x, envir = globalenv()))))[seq_len(num)]
}
