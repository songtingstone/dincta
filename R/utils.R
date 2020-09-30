#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom dplyr %>%
#' @examples
#' x <- 5 %>% sum(10)
#' 
#' @usage lhs \%>\% rhs
#' @return return value of rhs function. 
NULL


onehot <- function(x) {
    x_ = as.character(x)
    x_u = unique(x_)
    if(length(x_u)==1){
        res = matrix(1., nrow=length(x_),ncol=1)
        colnames(res) = as.character(x_u)[1]
    }else{
        res = data.frame(x_) %>% 
            tibble::rowid_to_column("row_id") %>% 
            dplyr::mutate(dummy = 1) %>% 
            tidyr::spread(x_, .data$dummy, fill = 0) %>% 
            dplyr::select(-.data$row_id) %>% as.matrix
    }
    return(res)
}

unknown_indicator <- function(x, unknown_type="unknown"){
    return(c(data.frame(x) == unknown_type))
}

scaleData <- function(A, margin = 1, thresh = 10) {
    if (!"dgCMatrix" %in% class(A))
        A <- methods::as(A, "dgCMatrix")
    
    if (margin != 1) A <- t(A)
    
    res <- scaleRows_dgc(A@x, A@p, A@i, ncol(A), nrow(A), thresh)
    if (margin != 1) res <- t(res)
    row.names(res) <- row.names(A)
    colnames(res) <- colnames(A)
    return(res)
}

moe_correct_ridge <- function(dinctaObj) {
    dinctaObj$moe_correct_ridge_cpp()
}


cluster <- function(dinctaObj, alpha, frequency_update_Phi_C) {
    if (dinctaObj$ran_init == FALSE) {
        stop('before clustering, run init_cluster')
    }
    dinctaObj$cluster_cpp(alpha,frequency_update_Phi_C)
}

# warm.up.iter.dincta, min.iter.dincta, max.iter.dincta,  alphas, verbose)
dinctate <- function(dinctaObj,  min_iter_dincta=3, iter_dincta=10,  alphas=NULL,frequency.update.Phi_C=NULL,
                     refine.Phi_C = TRUE,
                     keep.known.cell_type = TRUE,
                     max.times.refine.Phi_C = 3, 
                     epsilon.cell_type.changed = 1e-2,
                     select.refine.fraction = 0.2,
                     epsilon.cells.type.changed_frequence = 1e-7,
                     verbose=TRUE) {
    if (iter_dincta < 1) {
        return(0)
    }
    if (is.null(alphas) | length(alphas)<iter_dincta) {
        stop(gettextf('alphas cannot be NULL and must have at least length of max iter harmornys2 %d, but see %d ',
                      iter_dincta,length(alphas)))
    }
    times.refine.Phi_C = 1 
    for (iter in seq_len(iter_dincta)) {
        if (verbose) {
            message(gettextf('Dincta %d/%d', iter, iter_dincta))        
        }
        
        
        # STEP 1: do clustering
        err_status <- cluster(dinctaObj, alphas[iter], frequency.update.Phi_C[iter])
        if (err_status == -1) {
            stop('terminated by user')
        } else if (err_status != 0) {
            stop(gettextf('Dincta exited with non-zero exit status: %d', 
                          err_status))
        }
        
        # STEP 2: regress out covariates
        moe_correct_ridge(dinctaObj)
        
        # STEP 3: check for convergence
        if (dinctaObj$check_convergence(1)) {
            if(iter>=min_iter_dincta){
                if ( refine.Phi_C & 
                     (times.refine.Phi_C <= max.times.refine.Phi_C) & iter < iter_dincta ) {
                    if(dinctaObj$init_refine_Phi_C_and_check_u_convergence_cpp(
                        alphas[iter], keep.known.cell_type, epsilon.cell_type.changed,
                        select.refine.fraction,
                        epsilon.cells.type.changed_frequence)){
                        if (verbose) {
                            message(gettextf("Dincta refine Phi_C converged after %d times, and Dincta converged after %d iterations", 
                                             times.refine.Phi_C-1, iter)) 
                        }
                        return(0)
                    }else{
                        if (verbose) {
                            message(gettextf("Dincta refine Phi_C %d/%d", 
                                             times.refine.Phi_C, max.times.refine.Phi_C)) 
                        }
                        times.refine.Phi_C = times.refine.Phi_C + 1
                    } 
                }else{
                    if (verbose) {
                        message(gettextf("Dincta converged after %d iterations", 
                                         iter))    
                    }
                    return(0)
                }
                
            }else{
                if (verbose) {
                    message(gettextf("Dincta converged earlier after %d iterations,  continue run to the min iteration %d.", 
                                     iter,min_iter_dincta))    
                }
            }
            
        }
    }
}



dincta_classify <- function(dinctaObj, iter_dincta=10,  alphas=NULL,frequency.update.Phi_C=NULL,
                            refine.Phi_C = TRUE,
                            keep.known.cell_type = TRUE,
                            epsilon.cell_type.changed = 1e-2,
                            select.refine.fraction = 0.2,
                            epsilon.cells.type.changed_frequence = 1e-7,
                            verbose=TRUE) {
    if (iter_dincta < 1) {
        return(0)
    }
    if (is.null(alphas) | length(alphas)<iter_dincta) {
        stop(gettextf('alphas cannot be NULL and must have at least length of max iter harmornys2 %d, but see %d ',
                      iter_dincta,length(alphas)))
    }
    times.refine.Phi_C = 1 
    for (iter in seq_len(iter_dincta)) {
        if (verbose) {
            message(gettextf('Dincta %d/%d', iter, iter_dincta))        
        }
        
        
        # STEP 1: do clustering
        err_status <- cluster(dinctaObj, alphas[iter], frequency.update.Phi_C[iter])
        if (err_status == -1) {
            stop('terminated by user')
        } else if (err_status != 0) {
            stop(gettextf('Dincta exited with non-zero exit status: %d', 
                          err_status))
        }
        
        
        # STEP 2: check for convergence 
        # each number of dincta is the refine now 
        
        if ( refine.Phi_C & iter < iter_dincta ) {
            if(dinctaObj$init_refine_Phi_C_and_check_u_convergence_cpp(
                alphas[iter], keep.known.cell_type, epsilon.cell_type.changed,
                select.refine.fraction,
                epsilon.cells.type.changed_frequence)){
                if (verbose) {
                    message(gettextf("Dincta converged after %d iterations", 
                                     iter)) 
                }
                return(0)
            }else{
                times.refine.Phi_C = times.refine.Phi_C + 1
            } 
        }else{
            if (verbose) {
                message(gettextf("Dincta converged after %d iterations", 
                                 iter))    
            }
            return(0)
        }
        
        
        
        
    }
}


init_cluster <- function(dinctaObj, cluster_prior=NULL,alpha=1.) {
    if (dinctaObj$ran_setup == FALSE) {
        stop('before initializing cluster, run setup')
    }
    if (!is.null(cluster_prior)) {
        if (ncol(cluster_prior) != dinctaObj$N) {
            stop('cluster_prior must be defined by N cells')
        }
        if (nrow(cluster_prior) > dinctaObj$K) {
            stop('cluster_prior cannot contain more than K clusters')
        }
        C <- nrow(cluster_prior)
        dinctaObj$Y <- matrix(0, dinctaObj$d, dinctaObj$K)
        dinctaObj$Y[, seq_len(C)] <- compute_Y(dinctaObj$Z_cos,
                                               cluster_prior)
        dinctaObj$R <- matrix(0, dinctaObj$K, dinctaObj$N)
        dinctaObj$R[seq_len(nrow(cluster_prior)), ] <- cluster_prior
        
        
        ## if needed, initialize K-C clusters        
        if (C < dinctaObj$K) {
            Ynew <- t(stats::kmeans(t(dinctaObj$Z_cos), 
                                    centers = dinctaObj$K - C,
                                    iter.max = 25, nstart = 10)$centers)
            dinctaObj$Y[, seq(1+C, dinctaObj$K)] <- Ynew
        }
        
        dinctaObj$init_cluster_cpp(C,alpha)
    } else {
        dinctaObj$Y <- t(stats::kmeans(t(dinctaObj$Z_cos), 
                                       centers = dinctaObj$K,
                                       iter.max = 25, nstart = 10)$centers)
        dinctaObj$init_cluster_cpp(0,alpha)
    }
    
}


DinctaConvergencePlot <- function(
    dinctaObj, round_start=1, round_end=Inf, do_wrap=FALSE
) {  
    ## ignore initial value
    ## break down kmeans objective into rounds
    obj_fxn <- data.frame(
        kmeans_idx = Reduce(c, lapply(dinctaObj$kmeans_rounds, 
                                      function(rounds) {
                                          seq_len(rounds)
                                      })),
        dincta_idx = Reduce(c, lapply(
            seq_len(length(dinctaObj$kmeans_rounds)),
            function(i) {rep(i, dinctaObj$kmeans_rounds[i])})
        ),
        val = utils::tail(dinctaObj$objective_kmeans, -1)
    ) %>%
        dplyr::filter(.data$dincta_idx >= round_start) %>% 
        dplyr::filter(.data$dincta_idx <= round_end) %>% 
        tibble::rowid_to_column("idx") 
    
    
    plt <- obj_fxn %>% ggplot2::ggplot(ggplot2::aes(.data$idx, .data$val,
                                                    col = .data$dincta_idx)) + 
        ggplot2::geom_point(shape = 21) + 
        ggplot2::labs(y = "Objective Function", x = "Iteration Number")
    
    if (do_wrap) {
        plt <- plt + ggplot2::facet_grid(.~.data$dincta_idx, scales = 'free',
                                         space = 'free_x')
    } 
    return(plt)
}













