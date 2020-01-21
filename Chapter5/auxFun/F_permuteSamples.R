# Permute the samples for a null scenario
permuteSamples = function(mat){
        mat[sample(seq_len(nrow(mat))), ]
}