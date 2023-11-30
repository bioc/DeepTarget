## Function to calculate Drug Target Response (DTR).
DTR <- function(DN,GN,Pred,Exp,DRS,GES,CutOff=3,plot=TRUE ){
    if (nrow (Pred)==0) {stop("Pred should contain drug of interest")}
    L.c <- list(
        KO= colnames(GES),Drug_Prism= colnames(DRS),Expr= colnames(Exp))
    CL.M <- Reduce(intersect,L.c)
    if ( DN %in% Pred[,2]){
        if (plot){
            exp.g <- Exp[GN,CL.M]< CutOff
            ## True is low, False is high
            exp.g  <- factor(exp.g);
            colcode <- exp.g
            levels(colcode) <- c('High', 'Low')
            dat.c <- data.frame(
                GES.c = GES[GN,CL.M],
                DRS.c = DRS[Pred[DN,1],CL.M],
                Expr.c = colcode)
            ggplot(dat.c, 
                aes(x = dat.c[,1], y = dat.c[,2], color = factor(dat.c[,3])))+
                stat_smooth(method = 'lm')+
                theme_bw(base_size = 15)+
                labs(x = 'Viability after Cripspr Knockout', 
                    y = paste0('Response to ',Pred[DN,2]),
                    color =  paste0(GN,'\nExpression\nstatus'))+
                    stat_cor(label.y = c(1.3,1.4))}
    }else{
        stop("The drug of interest does not exist. Please double check")}
}

