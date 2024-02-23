# Function to obtain the datatype for Deeptarget.
Depmap2DeepTarget <- function(DataType="TPM", version="20Q4"){
    eh <- ExperimentHub()
    out <- query(eh, "depmap")
    ## check to see whether user enter correctly.
    Found.i <-  which ( out$title==paste0(DataType,"_",version))
    if (length(Found.i)==0){
        print ( "Only these datatypes and versions are available")
        print ( out$title[grep ( "mutationCalls",out$title)])
        print ( out$title[grep ( "TPM",out$title)])
        print ( out$title[grep ( "crispr",out$title)])
    }else{
        out.f <- out[[Found.i]]
        if (DataType=="TPM" ){
            out.d <- data.frame(ColumnName=out.f$depmap_id,RowName=out.f$gene_name,Value=out.f$rna_expression)
            out.m <- xtabs(Value ~ RowName + ColumnName, data = out.d)
        }
        if (DataType=="crispr") { 
            out.d <- data.frame(ColumnName=out.f$depmap_id,RowName=out.f$gene_name,Value=out.f$dependency)
            out.m <- xtabs(Value ~ RowName + ColumnName, data = out.d)}
        if (DataType=="mutationCalls"){
            mu_frame <- data.frame(ColumnName=out.f$depmap_id,
                                   RowName=out.f$gene_name)
            ## remove silent one.
            filter_var_class <- "Silent"
            f.i <- which(out.f$var_class %in% filter_var_class)
            mu_frame <- mu_frame[-f.i,]
            # Unique column and row names
            unique_cols <- unique(mu_frame$ColumnName)
            unique_rows <- unique(mu_frame$RowName)
            # Create an empty matrix
            out.m <- matrix(0, nrow = length(unique_rows), ncol = length(unique_cols),
                            dimnames = list(unique_rows, unique_cols))
            # Fill the matrix
            for(i in 1:nrow( out.m)) {
                out.m[mu_frame$RowName[i], mu_frame$ColumnName[i]] <- 1
            }
        }
        out.m
    }
}

