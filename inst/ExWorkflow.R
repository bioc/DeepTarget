## this is the example script if the data contains multiple assays for a drug.
rm(list=ls());
library ( DeepTarget)
data ("OntargetM")
drug.name <- c('ibrutinib','palbociclib')
dir.create ( 'Result')
## get the id,
S.Drug <- OntargetM$DrugMetadata$broad_id_trimmed [which (OntargetM$DrugMetadata$name %in% drug.name)]
sec.prism.f <- OntargetM$secondary_prism[which ( row.names(OntargetM$secondary_prism) %in% S.Drug), ]
KO.GES <- OntargetM$avana_CRISPR
## calculate the similarity between these assays with KO method.
List.sim <- NULL;
for ( i in 1:nrow(sec.prism.f)){
    DRS=as.data.frame(sec.prism.f[i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[i]
    out <- GetSim(row.names(sec.prism.f)[i],DRS=DRS, GES=KO.GES)
    List.sim [[length(List.sim) + 1]] <- out
}
names(List.sim) <- row.names(sec.prism.f)

####
metadata <- OntargetM$DrugMetadata
## get the similarity for known targeted gene from the drug ( if there are multiple targeted genes, get the most similarity)
DrugTargetSim <- PredTarget(Sim.GES.DRS=List.sim, D.M = metadata)
## get the similarity for the gene havnig the most similarty with the drug treatment.
DrugGeneMaxSim <- PredMaxSim(Sim.GES.DRS=List.sim, D.M = metadata)
## mutant interaction and expression interaction with the known targeted gene.
d.mt <- OntargetM$mutations_mat
d.expr <- OntargetM$expression_20Q4
out.MutantTarget <- NULL;
out.LowexpTarget <- NULL;
### if duplicated exists.
for ( i in 1:nrow(sec.prism.f)){
    identical (row.names(sec.prism.f)  , DrugTargetSim[,1])
    DRS=as.data.frame(sec.prism.f[i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[i]
    ## for mutant
    MutantInteract <- DoInteractMutant (Predtargets=DrugTargetSim[i,],Mutant=d.mt,DRS=DRS,GES=KO.GES)
    ## assign the estimate as strength, and P val from the interaction model.
    TargetMutSpecificity <- data.frame(MaxTgt_Inter_Mut_strength=sapply(MutantInteract, function(x) x[1]), MaxTgt_Inter_Mut_Pval=sapply(MutantInteract, function(x) x[2]))
    out.MutantTarget <- rbind (out.MutantTarget,TargetMutSpecificity)
    ## for expression.
    ExpInteract <- DoInteractExp (DrugTargetSim[i,],d.expr,DRS=DRS,GES=KO.GES,CutOff = 2)
    ## assign the estimate as strength, and P val from the interaction model.
    TargetExpSpecificity <- data.frame(MaxTgt_Inter_Exp_strength=sapply(ExpInteract, function(x) x[1]), MaxTgt_Inter_Exp_Pval=sapply(ExpInteract, function(x) x[2]))
    out.LowexpTarget <- rbind ( out.LowexpTarget,TargetExpSpecificity)
}

identical ( DrugTargetSim[,1],DrugGeneMaxSim[,1])

## whether interaction is true or false based on cut-off. estimate and p val from lm model
Whether_interaction_Ex_based= ifelse ( out.LowexpTarget$MaxTgt_Inter_Exp_strength <0 & out.LowexpTarget$MaxTgt_Inter_Exp_Pval <0.2,TRUE,FALSE)
## mutation interaction with P <0.1
predicted_resistance_mutation = ifelse ( out.MutantTarget$MaxTgt_Inter_Mut_Pval<0.1,TRUE,FALSE)
### if desired, save how many cellline has low expresion.
Pred.d <- cbind ( DrugTargetSim,DrugGeneMaxSim,out.MutantTarget,predicted_resistance_mutation, out.LowexpTarget,Whether_interaction_Ex_based)
Low.Exp = sapply(Pred.d[,3],function(x)errHandle(sum(d.expr[x,] < 2)) )
## save for later.
Pred.d$lowExpCount<-Low.Exp
##
### Obtain the similarity between the viablity scores from drug targeted treatment
## vs Gene effect score from KO method for the group that don't have primary targeted express.
## pred column 3.
## only use the rows with no NA.
idx <- which ( Pred.d$lowExpCount>0)
Pred.d.f <- Pred.d[idx ,]
Low.Exp.G = sapply(Pred.d.f[,3], function(x) errHandle(names(which(d.expr[x,]<2))))
identical ( names(Low.Exp.G),Pred.d.f[,3] )
## only perform the gene has at least some celllines having low exp
sim.LowExp <- NULL;
sec.prism.f.f <- sec.prism.f[idx,]
identical (row.names(sec.prism.f.f) ,Pred.d.f [,1])

for ( i in 1:nrow(Pred.d.f)){
    DRS.L= sec.prism.f.f[i,Low.Exp.G[[unlist(Pred.d.f[i,3])]]]
    DRS.L <- t(as.data.frame(DRS.L))
    row.names(DRS.L) <- Pred.d.f[i,1]
    out <- GetSim(Pred.d.f[i,1],DRS=DRS.L, GES=KO.GES)
    sim.LowExp [[length(sim.LowExp) + 1]] <- out
}
names(sim.LowExp) <-Pred.d.f[,1]
saveRDS(sim.LowExp,
        file = 'Result/similarity_KO_LowExp_DrugTreatment.RDS')

### plot and show the top 5 genes having the most corelation with drug when the primary is not expressed.
## make the function and use the sub function from dr. hu.

sim.LowExp.Strength=sapply(sim.LowExp, function(x) x[,2])
head(sim.LowExp.Strength)
sim.LowExp.Pval=sapply(sim.LowExp, function(x) x[,1])
head(sim.LowExp.Pval)
dim(sim.LowExp.Strength)
pdf ("Result/sim.low.exp.plot.pdf")
par(mar=c(4,4,5,2), xpd=TRUE, mfrow=c(2,2));
plotSim (dx=sim.LowExp.Pval,dy=sim.LowExp.Strength,clr=colorRampPalette(c("lightblue",'darkblue')), plot=TRUE)
dev.off();
## rcord these top 5 genes to the pred object.
L.topG <- NULL;
for ( i in 1:ncol(sim.LowExp.Strength)){
    top.5 <- names (sort(sim.LowExp.Strength[,i], decreasing=TRUE)[1:5])
    top.5  <- unlist(top.5);
    top.5    <- paste(top.5, collapse=" ");
    L.topG <- rbind ( L.topG, top.5 )
}
Pred.d$top5GeneWlowEx <- "NA"
Pred.d$top5GeneWlowEx [idx] <- L.topG
H.topG <- NULL;
simExp.Strength <- sapply(List.sim, function(x) x[,2])
dim(simExp.Strength)
for ( i in 1:ncol(simExp.Strength)){
    top.5 <- names (sort(simExp.Strength[,i], decreasing=TRUE)[1:5])
    top.5  <- unlist(top.5);
    top.5    <- paste(top.5, collapse=" ");
    H.topG <- rbind ( H.topG, top.5 )
}

Pred.d$top5Gene_Ex <- H.topG

##
write.csv (Pred.d,"./Result/Prediction_sim_KO_DrugTreatment.csv" )
## let plot for mutation =TRUE

## plot
DOI = 'dabrafenib'
GOI ='BRAF'
head(Pred.d)
## the first two cols is the drugs having two assays.
## we know that sec.prism.f has the same order with Pred.d.
which.mut <- which (Pred.d$predicted_resistance_mutation==TRUE);
## plot mutant.
## preaparing data
for ( i in 1:length(which.mut)){
    cr.i <- which.mut[i]
    DOI = Pred.d[cr.i,2]
    GOI =Pred.d[cr.i,3]
    DRS=as.data.frame(sec.prism.f[cr.i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[cr.i]
    head(DRS)
    pdf ( paste0("./Result/mutant_", row.names(Pred.d)[cr.i],".pdf"))
    out <- DMB (DN=DOI,GN=GOI,Pred=Pred.d[cr.i,],Mutant=d.mt,DRS= DRS,GES= KO.GES,plot=TRUE)
    print (out)
    dev.off();
}

## primary targeted.
which.exp <- which (Pred.d$Whether_interaction_Ex_based==FALSE);
for ( i in 1:length(which.exp )){
    cr.i <- which.exp[i]
    DOI = Pred.d[cr.i,2]
    GOI =Pred.d[cr.i,'MaxTargetName']
    DRS=as.data.frame(sec.prism.f[cr.i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[cr.i]
    head(DRS)
    pdf ( paste0("./Result/expr_", GOI,row.names(Pred.d)[cr.i],".pdf"))
    out <- DTR ( DN=DOI,GN=GOI,Pred=Pred.d[cr.i,], Exp=d.expr,DRS= DRS,GES=KO.GES,CutOff= 2)

    print (out)
    dev.off();
}

# based on the secondary targeted.
which.exp <- which (Pred.d$Whether_interaction_Ex_based==FALSE);
dim(Pred.d)
for ( i in 1:length(which.exp )){
    cr.i <- which.exp[i]
    DOI = Pred.d[cr.i,2]
    GOI =Pred.d[cr.i,'BestTargetGene']
    DRS=as.data.frame(sec.prism.f[cr.i,])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[cr.i]
    head(DRS)
    pdf ( paste0("./Result/expr_", GOI,row.names(Pred.d)[cr.i],".pdf"))
    out <- DTR ( DN=DOI,GN=GOI,Pred=Pred.d[cr.i,], Exp=d.expr,DRS= DRS,GES=KO.GES,CutOff= 2)
    print (out)
    dev.off();
}
### plot the correlation for the predited target.
DOI = 'atiprimod'
GOI = 'SOX10'
#GOI = 'MITF'
idx <- which( Pred.d[,2]==DOI)
dim(sec.prism.f)
identical ( row.names(sec.prism.f) , Pred.d$DrugID)
for ( i in 1:length(idx)){
    DRS=as.data.frame(sec.prism.f[idx[i],])
    DRS <- t(DRS)
    row.names(DRS) <- row.names(sec.prism.f)[idx[i]]
    head(DRS)
    pdf ( paste0("./Result/Cor_plot_Predicted_Of", GOI,'of_',row.names(Pred.d[idx[i],]),".pdf"));
    plotCor(DN=DOI,GN=GOI,Pred=Pred.d[idx[i],],DRS= DRS,GES= KO.GES,plot=TRUE);
    dev.off()
}



