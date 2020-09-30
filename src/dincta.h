#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <progress.hpp>


using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

using namespace std;
#include "dincta_types.h"
#include "utils.h"

class dincta_full;
RCPP_EXPOSED_CLASS(dincta_full)
  
  class dincta_full { 
  public: 
    /* CONSTRUCTORS etc */
    dincta_full(int __K);   
    
    void setup(MATTYPE& __Z, MATTYPE& __Phi, MATTYPE& __Phi_B, MATTYPE& __Phi_C,
               MATTYPE& __Phi_moe, ROWVECTYPE __u,
               float __sigma_entropy, 
               float __sigma_cell_type, 
               VECTYPE __theta_batch, 
               float __mu, 
               int __max_iter_kmeans,
               int __min_iter_cell_type,
               float __epsilon_kmeans,
               float __epsilon_dincta, 
               int __K, 
               int __C_known,
               int __C_residual,
               float __k_cluster_n_cells_outer_threshold,
               float __k_cluster_n_cells_inner_threshold,
               float __new_cell_type_prob_threshold,  
               float __new_cell_type_main_fraction,  
               float __cell_type_sample_fraction, 
               float __cell_type_min_cells,
               float __cell_type_eigval_threshold,
               float __new_cell_type_min_cells,
               float __centroid_cor_threshold,
               float __tau, float __block_size, 
               MATTYPE __lambda, bool __verbose, bool __print_loss);
    
    /* METHODS */
    void moe_correct_ridge_cpp();
    void init_cluster_cpp(unsigned Ch, float alpha);
    bool init_refine_Phi_C_and_check_u_convergence_cpp(
        float alpha, 
        bool keep_original_cell_type,
        float epsilon_cell_type_changed,
        float select_refine_fraction,
        float epsilon_changed_frequence);
    int cluster_cpp(float alpha, int frequency_update_Phi_C);
    
    void allocate_buffers();
    void compute_objective(float alpha); 
    int update_R(float alpha);
    void update_Phi_C(float alpha);
    bool check_convergence(int type);
    
    /* FIELDS */
    MATTYPE R, Z_orig, Z_corr, Z_cos, Y, Y_unnormed, Phi, Phi_B, Phi_C, Phi_moe;
    float sigma, sigma_entropy, sigma_cell_type, sigma_prior, tau;
    VECTYPE N_b, N_c, N_cb, Pr_b_given_c, theta_batch, theta_batch_vector;
    MATTYPE lambda; // diagonal MATTYPE matrix of ridge regression penalties
    vector<float> objective_dincta;
    vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, 
    objective_kmeans_cross, objective_kmeans_kl_cell_type;
    vector<float> centroid_max_cors,cell_type_min_cells_vector;
    vector<int> kmeans_rounds, update_Phi_C_rounds; 
    VECTYPE new_cell_type_rounds;
    
    
    float block_size, epsilon_kmeans, epsilon_dincta, merge_thresh_global, mu;
    int N, K, B, C, d, max_iter_kmeans, min_iter_cell_type, window_size; 
    // C = C_known+C_new+1+C_residual; maybe should keep the original
    // C_current = C_known+C_new+1; // will use to do computation for E, O
    // C_current_B = C_current * B; 
    int C_known, C_new, C_residual, C_current, C_current_B; 
    float  k_cluster_n_cells_outer_threshold, k_cluster_n_cells_inner_threshold,
    cell_type_sample_fraction, new_cell_type_main_fraction, 
    new_cell_type_prob_threshold, cell_type_eigval_threshold,cell_type_min_cells, 
    new_cell_type_min_cells, centroid_cor_threshold;
    
    // buffers
    // note O_c = E_c is N Pr(k,c)
    MATTYPE _scale_dist, dist_mat, O, O_b, E, E_c, E_t, dir_prior, Phi_Rk, theta_batch_matrix, A_U; 
    uvec update_order, cells_update, c_list, cells_type_unknown, cells_type_known,
    cells_type_known_original, cells_type_unknown_original;
    MATTYPE W, M, Phi_C_unknown, U, Ug, P; 
    ROWVECTYPE u_original, u_prior, u, u_scale;
    
    // flags
    bool ran_setup, ran_init, verbose, update_Phi_C_flag, Phi_C_changed,refine_Phi_C_flag, print_loss; 
    
  };


class dincta_marginal;
RCPP_EXPOSED_CLASS(dincta_marginal)
  
  class dincta_marginal { 
  public: 
    /* CONSTRUCTORS etc */
    dincta_marginal(int __K);    
    void setup(MATTYPE& __Z, MATTYPE& __Phi, MATTYPE& __Phi_B, MATTYPE& __Phi_C,
               MATTYPE& __Phi_moe, ROWVECTYPE __u,
               float __sigma_entropy, 
               float __sigma_cell_type, 
               VECTYPE __theta_batch, 
               float __mu, 
               int __max_iter_kmeans,
               int __min_iter_cell_type,
               float __epsilon_kmeans,
               float __epsilon_dincta, 
               int __K, 
               int __C_known,
               int __C_residual,
               float __k_cluster_n_cells_outer_threshold,
               float __k_cluster_n_cells_inner_threshold,
               float __new_cell_type_prob_threshold,  
               float __new_cell_type_main_fraction,  
               float __cell_type_sample_fraction, 
               float __cell_type_min_cells,
               float __cell_type_eigval_threshold,
               float __new_cell_type_min_cells,
               float __centroid_cor_threshold,
               float __tau, float __block_size, 
               MATTYPE __lambda, bool __verbose, bool __print_loss);
    
    
    /* METHODS */
    void moe_correct_ridge_cpp();
    void init_cluster_cpp(unsigned Ch, float alpha);
    bool init_refine_Phi_C_and_check_u_convergence_cpp(
        float alpha, 
        bool keep_original_cell_type,
        float epsilon_cell_type_changed,
        float select_refine_fraction,
        float epsilon_changed_frequence);
    int cluster_cpp(float alpha, int frequency_update_Phi_C);
    
    void allocate_buffers();
    void compute_objective(float alpha); 
    int update_R(float alpha);
    void update_Phi_C(float alpha);
    bool check_convergence(int type);
    
    /* FIELDS */
    MATTYPE R, Z_orig, Z_corr, Z_cos, Y, Y_unnormed, Phi, Phi_B, Phi_C, Phi_moe;
    float sigma, sigma_entropy, sigma_cell_type, sigma_prior, tau;
    VECTYPE N_b, N_c, N_cb, Pr_b_given_c, theta_batch, theta_batch_vector;
    MATTYPE lambda; // diagonal MATTYPE matrix of ridge regression penalties
    vector<float> objective_dincta;
    vector<float> objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy, 
    objective_kmeans_cross, objective_kmeans_kl_cell_type;
    vector<float> centroid_max_cors,cell_type_min_cells_vector;
    vector<int> kmeans_rounds, update_Phi_C_rounds; 
    VECTYPE new_cell_type_rounds;
    
    float block_size, epsilon_kmeans, epsilon_dincta, merge_thresh_global, mu;
    int N, K, B, C, d, max_iter_kmeans,min_iter_cell_type,  window_size; 
    // C = C_known+C_new+1+C_residual; maybe should keep the original
    // C_current = C_known+C_new+1; // will use to do computation for E, O
    // C_current_B = C_current * B; 
    int C_known, C_new, C_residual, C_current, C_current_B; 
    float   k_cluster_n_cells_outer_threshold, k_cluster_n_cells_inner_threshold,
    cell_type_sample_fraction, new_cell_type_main_fraction, new_cell_type_prob_threshold,
    cell_type_eigval_threshold,cell_type_min_cells, new_cell_type_min_cells, centroid_cor_threshold;
    
    // buffers
    // note O_c = E_c is N Pr(k,c)
    MATTYPE _scale_dist, dist_mat, O, O_b, E, E_b, E_c, E_t, dir_prior,
    Phi_Rk, theta_batch_matrix, A_U; 
    uvec update_order, cells_update, c_list, cells_type_unknown, cells_type_known,
    cells_type_known_original, cells_type_unknown_original;
    MATTYPE W, M, Phi_C_unknown, U, Ug, P; 
    ROWVECTYPE u_original, u_prior, u, u_scale; 
    
    // flags
    bool ran_setup, ran_init, verbose, update_Phi_C_flag, Phi_C_changed,refine_Phi_C_flag, print_loss;
  };


class dincta_classification;
RCPP_EXPOSED_CLASS(dincta_classification)
  
  class dincta_classification { 
  public: 
    /* CONSTRUCTORS etc */
    dincta_classification(int __K);    
    void setup(MATTYPE& __Z, MATTYPE& __Phi_C,
               ROWVECTYPE __u,
               float __sigma_entropy, 
               float __sigma_cell_type, 
               float __mu, 
               int __max_iter_kmeans,
               // int __min_iter_cell_type,
               float __epsilon_kmeans,
               float __epsilon_dincta,
               int __K, 
               int __C_known,
               int __C_residual,
               float __k_cluster_n_cells_outer_threshold,
               float __k_cluster_n_cells_inner_threshold,
               float __new_cell_type_prob_threshold,  
               float __new_cell_type_main_fraction,  
               float __cell_type_sample_fraction, 
               float __cell_type_min_cells,
               float __cell_type_eigval_threshold,
               float __new_cell_type_min_cells,
               float __centroid_cor_threshold,
               float __block_size, 
               bool __verbose, bool __print_loss);
    
    
    /* METHODS */
    void init_cluster_cpp(unsigned Ch, float alpha);
    bool init_refine_Phi_C_and_check_u_convergence_cpp(
        float alpha, 
        bool keep_original_cell_type,
        float epsilon_cell_type_changed,
        float select_refine_fraction,
        float epsilon_changed_frequence);
    int cluster_cpp(float alpha, int frequency_update_Phi_C);
    
    void allocate_buffers();
    void compute_objective(float alpha); 
    int update_R(float alpha);
    void update_Phi_C(float alpha);
    bool check_convergence(int type);
    
    /* FIELDS */
    MATTYPE R, Z_orig, Z_corr, Z_cos, Y, Y_unnormed, Phi_C;
    float sigma, sigma_entropy, sigma_cell_type, sigma_prior;
    VECTYPE N_c;
    vector<float> objective_dincta, objective_kmeans, objective_kmeans_dist, objective_kmeans_entropy,  objective_kmeans_kl_cell_type;
    vector<float> centroid_max_cors,cell_type_min_cells_vector;
    vector<int> kmeans_rounds, update_Phi_C_rounds; // OLD: Kb
    VECTYPE new_cell_type_rounds;
    
    float block_size, epsilon_kmeans, epsilon_dincta, mu;
    int N, K, C, d, max_iter_kmeans, min_iter_cell_type, window_size; 
    
    // C = C_known+1 + C_new +C_residual; maybe should keep the original
    // C_current = C_known+1+C_new; // will use to do computation for E, O
    // C_current_B = C_current * B; 
    int C_known, C_new, C_residual, C_current, C_current_B; 
    float   k_cluster_n_cells_outer_threshold, k_cluster_n_cells_inner_threshold,
    cell_type_sample_fraction, new_cell_type_main_fraction, new_cell_type_prob_threshold, 
    cell_type_eigval_threshold, cell_type_min_cells, new_cell_type_min_cells,centroid_cor_threshold;
    
    // buffers
    // note O_c = E_c is N Pr(k,c)
    MATTYPE _scale_dist, dist_mat, O, A_U; 
    uvec update_order, cells_update, c_list, cells_type_unknown, cells_type_known,
    cells_type_known_original, cells_type_unknown_original;
    MATTYPE M, Phi_C_unknown, U, Ug, P; 
    ROWVECTYPE u_original, u_prior, u, u_scale;
    
    // flags
    bool ran_setup, ran_init, verbose, update_Phi_C_flag, Phi_C_changed,refine_Phi_C_flag, print_loss; 
  };

