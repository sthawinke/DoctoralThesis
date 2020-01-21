#' Pseudo-F
pseudoF = function(coord, clusters, Dim = 1:2){
    coord = coord[,Dim, drop=FALSE]
    N = nrow(coord)
    a = length(unique(clusters))

    withinDist = sum(unlist(tapply(seq_len(nrow(coord)), clusters, function(x){
        dist(coord[x,])^2/length(x)
    })))

    overalDist = sum(dist(coord)^2)/N

    Fratio = (overalDist-withinDist)/withinDist * (a-1)/(N-a)
    return(Fratio)
}