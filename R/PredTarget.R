## this is used to pull out the targeted gene by the drug that has the max corelation.
## need to check this.

PredTarget <- function(Sim.GES.DRS,D.M){
    ## Make sure the user use the same meta data for the list of Sim.GES.DRS
    Drug.Id.i <- match(names(Sim.GES.DRS), D.M[,1])
    D.M.f <- D.M[Drug.Id.i,]
    ## based on the name of the drug.
    Splt.targets <- str_split(as.character(D.M.f[,3]), ", ")
    corrMatPval <- sapply(Sim.GES.DRS, function(x) x[,1])
    corrMat <- sapply(Sim.GES.DRS, function(x) x[,2])
    corrMatFDR <- sapply(Sim.GES.DRS, function(x) x[,3])
    Stat.Drug <- sapply(1:length(Splt.targets),function(x){
        ret_cor <- corrMat[
        match(Splt.targets[[x]],
        rownames(corrMat)), x]
        names(ret_cor) <- errHandle(Splt.targets[[x]])
        ret_cor
        })
    ## if there is one drug.
    if ( nrow(D.M.f)==1){
        ## add errHandle function to avoid the eror of returing -Inf
        Target.Max.cor <- errHandle(max(Stat.Drug,na.rm=TRUE ))
        Target.Max.Name <- names(Stat.Drug)[which.max(Stat.Drug)]
        ## drug Id will be the same as the metadata.
        Drug.id <- names(Sim.GES.DRS)
        names(Target.Max.cor) <- as.character(D.M.f[,2])
    }else{
        names(Stat.Drug) <- as.character(D.M.f[,2])
        Target.Max.cor <- sapply(Stat.Drug, function(x) max(x, na.rm=TRUE))
        Target.Max.Name <- sapply(Stat.Drug,function(x) names(which.max(x)))
        Drug.id <- D.M.f[match(names(Target.Max.cor), D.M.f$name),1]
    }
    Target.Max.Name[sapply(Target.Max.Name, length)==0] <- NA
    ### obtain drug ID ( make sure that we map the correct one)
    Target.Pred <- data.frame(
        DrugID = Drug.id,
        drugName=names(Target.Max.cor),
        MaxTargetName=unlist(Target.Max.Name),
        Maxcorr=Target.Max.cor)
    ### pull out p val and FDR.
    Target.Pred$KnownTargetCorrP <- sapply(1:nrow(Target.Pred), function(x)
        errHandle(corrMatPval[Target.Pred[x,3], Target.Pred[x,1]]) )
    # Known Target Significance - FDR corrected
    Target.Pred$KnownTargetCorrFDR <- sapply(1:nrow(Target.Pred), function(x)
        errHandle(corrMatFDR[Target.Pred[x,3], Target.Pred[x,1]]) )
    Target.Pred
}


