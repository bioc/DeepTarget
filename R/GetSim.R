GetSim <- function(DrugName,DRS,GES){
    ## make sure that users already use the common cell lines in both datasets
    common.c <- intersect(colnames(DRS),colnames(GES))
    DRS.c<- DRS[DrugName,common.c ]
    GES.c<- GES[,common.c]
    out.Sim.GES.DRS <- sapply(
        1:nrow(GES.c),
        function(x) unlist(cor.test_trimmed_v0.default(DRS.c, GES.c[x,])))
    out.Sim.GES.DRS <- t(out.Sim.GES.DRS)
    row.names(out.Sim.GES.DRS  ) <- row.names(GES.c)
    FDR <- fdrCor(out.Sim.GES.DRS[,1])
    out.Sim.GES.DRS <- cbind(out.Sim.GES.DRS, FDR=FDR)
    out.Sim.GES.DRS
}
