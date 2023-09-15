## perform the interaction based on cut-off
DoInteractMutant <- function(Predtargets,Mutant,DRS,GES){
    #getting only the overlapped samples
    L.c <- list( KO= colnames(GES),
                Drug_Prism= colnames(DRS),
                Mutation= colnames(Mutant))
    CL.M <- Reduce(intersect,L.c)
    interactFeatures <- lapply(1:nrow(Predtargets), function(x){
            ## drug response based on the drug ID
            DRS.f <- errHandle(DRS[Predtargets[x,1],CL.M])
            ## Gene effect scores based on the targeted gene
            GES.f <- errHandle(GES[Predtargets[x,3],CL.M])
            Mutant.Group <- errHandle(Mutant[Predtargets[x,3],CL.M] )
            ## True is low, False is high.
            errHandle(summary(
                lm(DRS.f ~ GES.f * Mutant.Group))$coefficients[4,c(1,4)])
            })
    names(interactFeatures) <- Predtargets$drugName
    interactFeatures
}


