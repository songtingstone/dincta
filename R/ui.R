#' Main Dincta interface
#' 
#' Use this to run the Dincta algorithm on gene expression or PCA matrix. 
#' 
#' @param data_mat Matrix of normalized gene expession (default) or PCA 
#' embeddings (see do_pca). 
#' Cells can be rows or columns. 
#' @param meta_data Dataframe with variables to integrate.
#' @param vars_use The variable(s) in meta_data 
#' will be removed  batch effects(character vector).
#' @param cell_type_use The variable 
#' contains cell type information (character vector) in meta_data.
#' @param R.cross.entropy.type The cross entropy type of clustering loss. Must be 'b' or "bc".
#' If it is 'b', using the O_kb/E_kb in the cross entropy loss. If it is 'bc', using the O_kcb/E_kcb 
#' in the cross entropy loss. Default R.cross.entropy.type = 'bc'.
#' @param cell_type_unknown The name of the unknown cell type (character), 
#' Default cell_type_unknow = "unknown".
#' @param Phi_C The initial cell type assignment matrix, Optional, can be NULL. If provided, then its dim(Phi_C) ==c(C,N).
#' Will update the cell type assignment Phi_C[,i] of cell i with meta_data[[cell_type_use]][i] == cell_type_unknown;
#' and keep the cell type assignments Phi_C[,i] of cells with known cell type,
#' i.e. meta_data[[cell_type_use]][i] != cell_type_unknown.
#' @param do_pca Whether to perform PCA on input matrix. 
#' @param npcs If doing PCA on input matrix, number of PCs to compute. 
#' @param theta.batch Diversity clustering penalty parameter. Specify for each
#'  variable in vars_use Default theta.batch=2. theta.batch=0 does not encourage any 
#'  diversity. Larger values of theta result in more diverse clusters. 
#' @param alphas The scaling factor for unknown cell type cells in each round of Dincta. (float or float vector). 
#' Between 0 to 1, default 0.5. 
#' @param lambda Ridge regression penalty parameter. Specify for each variable
#'  in vars_use. 
#' Default lambda=1. Lambda must be strictly positive. Smaller values result 
#' in more aggressive correction. 
#' @param sigma.entropy  Entropy constriant penalty parameter. Default sigma.entropy=0.1.
#' @param sigma.cell_type Cell type constraint penalty parameter. Default sigma.cell_type=0.1. 
#'  Larger values of sigma.cell_type result in more strict cell type constraints. 
#'  Suggestion sigma.cell_type <= sigma.entropy.
#' @param mu A scalar to control the sparsity in R update, defualt mu = 0.01
#' @param nclust Number of clusters in model. nclust=1 equivalent to simple 
#' linear regression. 
#' @param n.cell_type.residual The max number of new cell types occurring in the data allowed. 
#' @param k_cluster.n.cells.outer.threshold The threshold (scalar) which determines the k cluster is the new cluster or not.
#' when [R(mathcal I) mathbf 1]_k < k_cluster.n.cells.outer.threshold will be regarded as the new cluster.
#' @param k_cluster.n.cells.inner.threshold The threshold (scalar) which determines  whether select k cluster to compute the A or not.
#' when [R mathbf 1]_k > k_cluster.n.cells.inner.threshold will be selected.
#' @param new.cell_type.prob.threshold The threshold which determine whether the cells belong to the new cell type or not.
#' when [sum_{k in mathcal K_{new}} R_ki < new.cell_type.prob.threshold], then the cell i will be regard as the candidate 
#' cells belong to the new cell type.
#' @param cell_type.sample.fraction The fraction of the all cells will be selected as the sample set.
#' @param new.cell_type.main.fraction The fraction of the candidate cells belong to the new cell type will be selected as the main set.
#' @param  cell_type.min.cells the min number of  cells that  a cell type should have. If one candiate cell type determined by [lambda(A)(k)]  have the number
#' the effective cells in the sample set less than cell_type.min.cells, then it will merge to its nearby cell type. 
#' @param cell_type.eigval.threshold The threshold determine the number of the cell types. 
#' When [lambda(A)(k) < cell_type.eigval.threshold] then $k$ will be count as one cell type.
#' @param  new.cell_type.min.cells the min number of  cells that  a new cell type should have. If one new cell type determined by [lambda(A)(k)]  have the number
#' the effective cells in the sample set less than new.cell_type.min.cells, then it will merge to its nearby cell type. 
#' @param centroid.cor.threshold A threshold to determine whether merge two cell type. 
#' Suppose that Y_c1 is the centroid of the cell type c1,   Y_c2 is the centroid of the cell type c2, if <Y_c1/||Y_c1||,Y_c2/||Y_c2/||> < centroid.cor.threshold; 
#' then c1 amd c2 will be merged to one cell type. 
#' @param tau Protection against overclustering small datasets with large ones.
#'  tau is the expected number of cells per cluster. 
#' @param block.size What proportion of cells to update during clustering.
#'  Between 0 to 1, default 0.05. Larger values may be faster but less accurate
#' @param min.iter.dincta Maximum number of rounds of dincta to prevent early convegent of R.  Default min.iter.dincta = 5.
#' @param max.iter.dincta Maximum number of rounds to run Dincta. One round
#'  of Dincta involves one clustering and one correction step. Default max.iter.dincta = 10. 
#' @param max.iter.cluster Maximum number of rounds to run clustering at each 
#' round of Dincta. Default  max.iter.cluster = 200.
#' @param frequency.update.Phi_C Upating Phi_C per frequency.update.Phi_C rounds at clustering.
#' @param refine.Phi_C Whether to print refine cell type assignments Phi_C. TRUE to refine
#' @param keep.known.cell_type Whether to keep the original known cell type assignment unchanged when refining Phi_C.
#' @param max.times.refine.Phi_C The max times to refine Phi_C.It will run dincta iteration until convergent in each time.
#' @param select.refine.fraction The fraction of all cells to refine the cell type assignment Phi_C. 
#' @param epsilon.cluster Convergence tolerance for clustering round of 
#' Dincta. Set to -Inf to never stop early. Default epsilon.cluster = 1e-5.
#' @param epsilon.dincta Convergence tolerance for Dincta. Set to -Inf to
#'  never stop early. Default epsilon.dincta = 1e-4.
#' @param epsilon.cell_type.changed The threshold to identify the cell type changed or not. Specifically, 
#' let Psi_C = Phi_C R^T diag^{-1} (R 1_N) R  be the cell type assignment compute from cluster assignments R, 
#' if sum(Phi_C[,i] log (Phi_C[,i]/Psi_C[,i]))  > epsilon.cell_type.changed, then cell i will be regarded  as  the cell with unknown cell type. 
#' @param epsilon.cells.type.changed_frequence  Convergence tolerance for refining  cell type assignment matrix Phi_C, 
#' when the fraction of the number of cells  which have the unknown cell types less 
#' than the epsilon.cells.type.changed_frequence, stop refining cell type assignment matrix Phi_C.      . 
#' @param plot_convergence Whether to print the convergence plot of the 
#' clustering objective function. TRUE to plot, FALSE to suppress. This can be
#'  useful for debugging. 
#' @param return_object (Advanced Usage) Whether to return the Dincta object 
#' or only the corrected PCA embeddings. 
#' @param reference_values (Advanced Usage) Defines reference dataset(s). 
#' Cells that have batch variables values matching reference_values will not 
#' be moved.  
#' @param cluster_prior (Advanced Usage) Provides user defined clusters for 
#' cluster initialization. If the number of provided clusters C is less than K, 
#' Dincta will initialize K-C clusters with kmeans. C cannot exceed K.  
#' @param verbose Whether to print progress messages. TRUE to print, 
#' FALSE to suppress.
#' @param print.loss Whether to print loss in each step. TRUE to print, 
#' FALSE to suppress.
#'
#' @return By default, a list consist of the matrix with corrected PCA embeddings, meta_data, and the Phi_C.
#' The meta_data  contains the predict cell type of the unknown cell type 
#' and the prob  that the cell belong to the predict cell type.
#' The Phi_C is the cell type assign matrix
#'   If return_object 
#' is TRUE, returns a list consist of the full Dincta object (R6 reference class type) , meta_data and the Phi_C.
#' The meta_data  and Phi_C are same as described in the default case. 
#' 
#' @export 
#' 
#' @examples
#' 
#' 
#' ## By default, Dincta inputs a normalized gene expression matrix
#' \dontrun{
#' dincta_res <- DinctaMatrix(exprs_matrix, meta_data, 'dataset','cell_type')
#' dincta_embedding <- dincta_res$dincta_embedding
#' meta_data <- dincta_res$meta_data
#' Phi_C <- dincta_res$Phi_C
#' }
#' 
#' ## Dincta can also take a PCA embeddings matrix
#' data(cell_lines)
#' pca_matrix <- cell_lines$scaled_pcs
#' meta_data <- cell_lines$meta_data
#' dincta_res <- DinctaMatrix(pca_matrix, meta_data, 'dataset', 'cell_type', 
#'                                     do_pca=FALSE)
#' ## dincta_res <- DinctaMatrix(pca_matrix, meta_data, 'dataset', 'partial_unknown_cell_type',  do_pca=FALSE)
#'                                    
#' ## dincta_res <- DinctaMatrix(pca_matrix, meta_data, 'dataset', 'whole_unknown_cell_type', do_pca=FALSE)
#'                                     
#' dincta_embedding <- dincta_res$dincta_embedding
#' meta_data <- dincta_res$meta_data
#' Phi_C <- dincta_res$Phi_C                                    
#' 
#' ## Output is a matrix of corrected PC embedding and predict cell type in the meta data
#' dim(dincta_embedding)
#' dincta_embedding[seq_len(5), seq_len(5)]
#' head(meta_data)
#' Phi_C[, seq_len(5)]
#' 
#' ## Finally, we can return an object with all the underlying data structures and
#' ## also the  predict cell type in the meta data
#' dincta_res <- DinctaMatrix(pca_matrix, meta_data, 'dataset', 'cell_type', 
#'                                     do_pca=FALSE, return_object=TRUE)
#' ## dincta_res <- DinctaMatrix(pca_matrix, meta_data, 'dataset', 'partial_unknown_cell_type',  do_pca=FALSE, return_object=TRUE)
#' ## dincta_res <- DinctaMatrix(pca_matrix, meta_data, 'dataset', 'whole_unknown_cell_type',  do_pca=FALSE, return_object=TRUE)                                
#'                                     
#' dincta_object <- dincta_res$dincta_object
#' meta_data <- dincta_res$meta_data
#' Phi_C <- dincta_res$Phi_C                                    
#' dim(dincta_object$Y) ## cluster centroids
#' dim(dincta_object$R) ## soft cluster assignment
#' dim(dincta_object$Z_corr) ## corrected PCA embeddings
#' head(meta_data) ## meta data contains the predict cell type of the unknown cell type  and the prob  that the cell belong to the predict cell type.
#' Phi_C[, seq_len(5)] ## cell type assign matrix
#' 
DinctaMatrix <- function(
    data_mat, meta_data, vars_use, cell_type_use=NULL, 
    R.cross.entropy.type ='bc',
    cell_type_unknown="unknown", 
    Phi_C = NULL,
    do_pca = TRUE,
    npcs = 20, 
    theta.batch = NULL, 
    alphas = 0.5,
    lambda = NULL, sigma.entropy = 0.1, sigma.cell_type = 0.1,
    mu =0.01, 
    nclust = NULL, 
    n.cell_type.residual=16,
    k_cluster.n.cells.outer.threshold=6,
    k_cluster.n.cells.inner.threshold=1,
    new.cell_type.prob.threshold=0.8,
    cell_type.sample.fraction=0.9,
    new.cell_type.main.fraction = 0.9,
    cell_type.eigval.threshold=0.9,
    cell_type.min.cells=10,
    new.cell_type.min.cells = 6,
    centroid.cor.threshold= 0.6,
    tau = 0, block.size = 0.05,
    min.iter.dincta = 3,
    min.iter.cell_type =3, 
    max.iter.dincta = 20, 
    max.iter.cluster = 200,
    frequency.update.Phi_C =NULL,
    refine.Phi_C = FALSE,
    keep.known.cell_type = TRUE,
    max.times.refine.Phi_C = 3,
    select.refine.fraction = 0.2,
    epsilon.cluster = 1e-5,
    epsilon.dincta = 1e-4, 
    epsilon.cell_type.changed = 1e-2,
    epsilon.cells.type.changed_frequence = 1e-7,
    plot_convergence = FALSE, return_object = FALSE, 
    reference_values = NULL, cluster_prior = NULL, 
    verbose = TRUE, print.loss = FALSE
) {
    if (!(R.cross.entropy.type %in% c('b','bc'))){
        stop('R.cross.entropy.type  must be either "b" 
             (using Pr(k,b) to compute R cross entropy loss ) 
             or "bc" (using Pr(k,b,c) to compute R cross entropy loss)')
    }
    # the meta_data[[cell_type_use]] can only include the unknown cell type information
    # and provide the cell type informaton in Phi_C, and if provide Phi_C, it has a higher priority .
    if (!(is(meta_data, 'data.frame') | is(meta_data, 'DataFrame'))) {
        #     if (!c('data.frame', '') %in% class(meta_data)) {
        if (length(meta_data) %in% dim(data_mat)) {
            vars_use <- 'batch_variable'
            cell_type_use <- 'cell_type'
            cell_type = rep("unknown", length(meta_data))
            meta_data <- data.frame(batch_variable = meta_data, cell_type = cell_type)
        } else {
            stop('meta_data must be either a data.frame or a vector with batch 
                values for each cell')
        }
    }
    
    if (is.null(vars_use) | any(!vars_use %in% colnames(meta_data))) {
        msg <- gettextf('must provide variables names (e.g. vars_use=%s)', 
                        sQuote('stim'))
        stop(msg)
    }
    if (is.null(cell_type_use) | any(!cell_type_use %in% colnames(meta_data))) {
        msg <- gettextf('must provide cell_type names (e.g. cell_type_use=%s)', 
                        sQuote('cell_type'))
        stop(msg)
    }
    
    if (do_pca) {
        if (ncol(data_mat) != nrow(meta_data)) {
            data_mat <- Matrix::t(data_mat)
        }
        
        pca_res <- data_mat %>%
            scaleData() %>% 
            irlba::prcomp_irlba(n = npcs, retx = TRUE, center = FALSE, 
                                scale. = FALSE)
        data_mat <- pca_res$rotation %*% diag(pca_res$sdev)
    } 
    
    N <- nrow(meta_data)
    cells_as_cols <- TRUE
    if (ncol(data_mat) != N) {
        if (nrow(data_mat) == N) {
            data_mat <- t(data_mat)
            cells_as_cols <- FALSE
        } else {
            stop("number of labels do not correspond to number of 
                samples in data matrix")
        }
    }
    
    if (is.null(nclust)) {
        nclust <- min(round(N / 30), 100)
    }
    n_alphas = length(alphas)
    if (n_alphas==1){
        alphas = rep(alphas,max.iter.dincta)
    }else if (n_alphas < max.iter.dincta){
        alphas = c(alphas, rep(alphas[n_alphas],max.iter.dincta - min.iter.dincta ))
    }
    if (is.null(theta.batch)) {
        theta.batch<- rep(2, length(vars_use))
        # theta.batch<-2
    }
    if (is.null(lambda)) {
        lambda <- rep(1, length(vars_use))
    }    
    
    if (is.null(mu)) {
        lambda <- rep(0.05, length(vars_use))
    }
    
    ## Pre-compute some useful statistics
    Phi_B<- Reduce(rbind, lapply(vars_use, function(var_use) {
        t(onehot(meta_data[[var_use]]))
    }))
    
    B = dim(Phi_B)[1]
    # can provide Phi_C, e.g from DinctaClassfication results
    if(is.null(Phi_C)){
        Phi_C<- Reduce(rbind, lapply(cell_type_use, function(ct_use) {
            t(onehot(meta_data[[ct_use]]))
        })) 
    }
    
    # reorder the rows of Phi_C, make sure that "unknown" cell type in the last row
    if(length(rownames(Phi_C)) >1){
        sort_res = sort(rownames(Phi_C), index.return =TRUE )
        ix = sort_res$ix
        ix_reorder = c(ix[sort_res$x != cell_type_unknown], ix[sort_res$x == cell_type_unknown])
        Phi_C = Phi_C[ix_reorder,]
    }
    
    # unkown cell type indicator
    # current only support one cell type 
    u =  Reduce(c, lapply(cell_type_use, function(ct_use) {
        unknown_indicator(meta_data[[ct_use]],cell_type_unknown)
    }))
    if(length(u)!=N){
        stop("only support for one cell type used! ")
    }else if(sum(u)>0 & (! cell_type_unknown %in% rownames(Phi_C))){
        stop(paste("the unknown cell type '", cell_type_unknown, 
                   "' not in the rownames of Phi_C, whose rownames is: ", paste(rownames(Phi_C),collapse=", "), sep = "" ))
    }
    
    
    C = dim(Phi_C)[1]
    if(sum(u)==0){
        n.cell_type.known = C
    }else{
        n.cell_type.known = C-1
    }
    
    # frequency.update.Phi_C 
    if(is.null(frequency.update.Phi_C)){
        frequency.update.Phi_C = rep(max.iter.cluster+1, max.iter.dincta) # only once in each iteration
        if(sum(u)<N*0.9){
            frequency.update.Phi_C[min.iter.dincta:max.iter.dincta] = 3
        }else{
            frequency.update.Phi_C[min.iter.dincta:max.iter.dincta] = 5
        } 
    }else if(length(frequency.update.Phi_C)==1){
        tmp = frequency.update.Phi_C
        frequency.update.Phi_C = rep(max.iter.cluster+1, max.iter.dincta) # only once in each iteration
        frequency.update.Phi_C[min.iter.dincta:max.iter.dincta] = tmp 
    }else if(length(frequency.update.Phi_C) < max.iter.dincta){
        tmp =  frequency.update.Phi_C
        n_tmp = length(tmp)
        frequency.update.Phi_C = c(tmp, rep(tmp[n_tmp], max.iter.dincta - n_tmp ))
    }
   
    Phi<- Reduce(rbind, lapply(seq_len(C), function(i) {
        matrix(rep(Phi_C[i,],B),nrow=B,ncol=length(Phi_C[i,]),byrow = T)*Phi_B
    }))
    
    B_vec <- Reduce(c, lapply(vars_use, function(var_use) {
        length(unique(meta_data[[var_use]]))
    }))
    
    # simullate annealing, length B
    theta.batch <- Reduce(c, lapply(seq_len(length(B_vec)), function(b) 
        rep(theta.batch[b], B_vec[b])))
    
    lambda <- Reduce(c, lapply(seq_len(length(B_vec)), function(b) 
        rep(lambda[b], B_vec[b])))
    # using 1e-8  in the intercept term to avoid R_k is zeros
    lambda_mat <- diag(c(1e-8, lambda))
    
    ## TODO: check that each ref val matches exactly one covariate
    ## TODO: check that you haven't marked all cells as reference! 
    if (!is.null(reference_values)) {
        idx <- which(row.names(Phi_B) %in% reference_values)
        cells_ref <- which(colSums(Phi_B[idx, , drop = FALSE] == 1) >= 1)
        b_keep <- which(!row.names(Phi_B) %in% reference_values)
        Phi_moe <- Phi_B[b_keep, , drop = FALSE]
        Phi_moe[, cells_ref] <- 0
        
        Phi_moe <- rbind(rep(1, N), Phi_moe)
        lambda_mat <- lambda_mat[c(1, b_keep + 1), c(1, b_keep + 1)]    
    } else {
        Phi_moe <- rbind(rep(1, N), Phi_B)
    }
    
    ## RUN Dincta
    if(R.cross.entropy.type == 'b'){
        dinctaObj <- new(dincta_marginal, 0) ## 0 is a dummy variable - will change later
    }else{
        dinctaObj <- new(dincta_full, 0)
    }
    dinctaObj$setup(
        data_mat, Phi, Phi_B, Phi_C, Phi_moe, 
        u, 
        sigma.entropy, 
        sigma.cell_type,
        theta.batch,
        mu,
        max.iter.cluster,
        min.iter.cell_type,
        epsilon.cluster,
        epsilon.dincta,
        nclust, 
        n.cell_type.known,
        n.cell_type.residual,
        k_cluster.n.cells.outer.threshold,
        k_cluster.n.cells.inner.threshold,
        new.cell_type.prob.threshold, 
        new.cell_type.main.fraction,
        cell_type.sample.fraction,
        cell_type.min.cells,
        cell_type.eigval.threshold,
        new.cell_type.min.cells,
        centroid.cor.threshold,
        tau, block.size, lambda_mat, verbose, print.loss
    )
    if((sum(dinctaObj$u_original) ==0) & keep.known.cell_type ){
        refine.Phi_C = FALSE
    }
    init_cluster(dinctaObj, cluster_prior,alphas[1])
    
    dinctate(dinctaObj, min.iter.dincta, max.iter.dincta, alphas,frequency.update.Phi_C,
             refine.Phi_C,
             keep.known.cell_type,
             max.times.refine.Phi_C,
             epsilon.cell_type.changed,
             select.refine.fraction,
             epsilon.cells.type.changed_frequence,
             verbose)
    if (plot_convergence) graphics::plot(DinctaConvergencePlot(dinctaObj))
    
    
    # get the predict cell type
    cell_types_original = rownames(Phi_C)
    n_cell_types_original = dim(Phi_C)[1]
    C_known =  dinctaObj$C_known
    C_new = dinctaObj$C_new
    n_cell_types_predict = n_cell_types_original + C_new
    Phi_C_predict = dinctaObj$Phi_C[1:n_cell_types_predict,]
    
    #n_cell_types_predict = dim(Phi_C_predict)[1]
    cell_types_predict = cell_types_original
    if(n_cell_types_original< n_cell_types_predict){
        new_cell_types_name = paste("new_cell_type_", 1:C_new, sep="")
        cell_types_predict = c(cell_types_predict, new_cell_types_name)
    }
    colnames(Phi_C_predict)= colnames(data_mat)
    rownames(Phi_C_predict)= cell_types_predict
    cell_type_predict_index = apply(Phi_C_predict, 2, which.max)
    cell_type_predict_prob = apply(Phi_C_predict,2, max)
    meta_data[['cell_type_predict']]  = as.character(lapply(cell_type_predict_index, function(cti){
        as.character(cell_types_predict[cti])
    }))
    meta_data[['cell_type_predict_prob']] = c(cell_type_predict_prob)
    if (return_object) {
        result = list()
        result$dincta_object = dinctaObj
        result$meta_data = meta_data
        result$Phi_C = Phi_C_predict
        return(result)
    } else {
        res <- as.matrix(dinctaObj$Z_corr)
        row.names(res) <- row.names(data_mat)
        colnames(res) <- colnames(data_mat)
        if (!cells_as_cols) 
            res <- t(res)
        result = list()
        result$dincta_embedding = res 
        result$meta_data = meta_data
        result$Phi_C = Phi_C_predict
        return(result)      
    }
}

