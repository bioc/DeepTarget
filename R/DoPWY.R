## Function to execute the pathway analysis for a set of targets.
DoPWY <- function(Sim.GES.DRS,D.M){
    PwyTargeted <- split(D.M[,3], D.M[,5])
    PwyTargeted.unl <- lapply(PwyTargeted, unlist)
    PwyTargeted.c <- lapply(PwyTargeted.unl, function(y){
        unique(unlist(lapply(as.character(y), function(x)
            unlist(strsplit(x, ', ')))))
    })
    PwyTargeted.c <- bplapply( PwyTargeted.c, na.omit)
    F.V <- nrow(Sim.GES.DRS[[1]])
    corrMat <- vapply(Sim.GES.DRS, function(x) x[,2],numeric(F.V))
    DoPWYEnr <- bplapply(seq_len(ncol(corrMat)), function(x) {
        corrMat.u <- unlist(corrMat[,x])
        names(corrMat.u) <- rownames(corrMat)
        corrMat.u.O <- sort(corrMat.u, decreasing = TRUE)
        gseaEnr <- fgsea(
            pathways = PwyTargeted.c,stats = corrMat.u.O,minSize=1,maxSize=100)
        gseaEnr})
    names(DoPWYEnr) <- colnames(corrMat)
    DoPWYEnr
}

