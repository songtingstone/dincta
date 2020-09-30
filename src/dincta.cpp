#include "dincta.h"



// NOTE: This is a dummy constructor, needed by Rcpp
dincta_full::dincta_full(int __K): K(__K) {}



void dincta_full::setup(MATTYPE& __Z, MATTYPE& __Phi, MATTYPE& __Phi_B, MATTYPE& __Phi_C,
                        MATTYPE& __Phi_moe, ROWVECTYPE __u,
                        float __sigma_entropy, 
                        float __sigma_cell_type,
                        VECTYPE __theta_batch, // vector with size B
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
                        MATTYPE __lambda, bool __verbose, bool __print_loss) {
  
  Z_corr = MATTYPE(__Z);
  Z_orig = MATTYPE(__Z);
  Z_cos = MATTYPE(Z_orig);
  cosine_normalize(Z_cos, 0, true); // normalize columns
  
  
  Phi_B = __Phi_B;
  Phi_moe = __Phi_moe;
  u = __u;
  u_prior = 16. - u; // just make the first check_u_convergence return false
  u_original = u; 
  N = Z_corr.n_cols;
  B = Phi_B.n_rows;
  d = Z_corr.n_rows; 
  window_size = 3;
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_dincta = __epsilon_dincta;
  
  lambda = __lambda;
  sigma_entropy = __sigma_entropy;
  sigma_cell_type = __sigma_cell_type;
  sigma =  sigma_entropy + sigma_cell_type;
  mu = __mu; 
  sigma_prior = sigma;
  block_size = __block_size;
  K = __K;
  C_known = __C_known;
  C_residual = __C_residual; // use for the new cell types
  k_cluster_n_cells_outer_threshold = __k_cluster_n_cells_outer_threshold;
  k_cluster_n_cells_inner_threshold = __k_cluster_n_cells_inner_threshold;
  new_cell_type_prob_threshold = __new_cell_type_prob_threshold;
  cell_type_sample_fraction = __cell_type_sample_fraction;
  new_cell_type_main_fraction = __new_cell_type_main_fraction;
  cell_type_eigval_threshold = __cell_type_eigval_threshold;
  cell_type_min_cells = __cell_type_min_cells;
  new_cell_type_min_cells = __new_cell_type_min_cells;
  centroid_cor_threshold = __centroid_cor_threshold;
  
  tau = __tau;
  // full cell type informaton 
  if(sum(u)==0){
    C_current = C_known;
    C_current_B = C_current * B;
    C_new = 0;
    C_residual = 0;
    C = C_known;
    Phi_C = __Phi_C;
    Phi = __Phi;
  }else{
    C_current = C_known + 1;
    C_current_B = C_current * B;
    C_new = 0;
    C = C_known + 1 + C_residual;
    Phi_C = zeros<MATTYPE>(C, N); 
    Phi_C.head_rows(C_current) = __Phi_C;
    Phi = zeros<MATTYPE>(C*B, N);
    Phi.head_rows(C_current*B) = __Phi;
  }
  max_iter_kmeans = __max_iter_kmeans;
  min_iter_cell_type = __min_iter_cell_type;
  verbose = __verbose;
  print_loss = __print_loss;
  
  theta_batch = __theta_batch;
  //theta_batch_vector = theta_batch*ones<VECTYPE>(C*B);
  theta_batch = repmat(theta_batch, C,1);// unnormalize, length C*B
  theta_batch_vector = theta_batch; // may be normalized by N_cb, length C*B
  //theta_batch_matrix = arma::repmat(theta_batch.t(), K, 1);
  
  // unknown cell type index 
  cells_type_unknown = find(u>0); 
  cells_type_unknown_original = cells_type_unknown;
  cells_type_known = find(u==0); 
  cells_type_known_original = cells_type_known;
  update_Phi_C_flag = cells_type_unknown.size()>0;
  
  allocate_buffers();
  ran_setup = true;
  refine_Phi_C_flag = false;
}


void dincta_full::allocate_buffers() {
  //_scale_dist = zeros<MATTYPE>(K, N);    
  dist_mat = zeros<MATTYPE>(K, N); 
  O = zeros<MATTYPE>(K, C*B);
  E = zeros<MATTYPE>(K, C*B); 
  E_t = zeros<MATTYPE>(K, C*B); 
  E_c = zeros<MATTYPE>(K,C);
  W = zeros<MATTYPE>(B + 1, d); 
  Phi_Rk = zeros<MATTYPE>(B + 1, N);
  M =  zeros<MATTYPE>(K,C);
  U = zeros<MATTYPE>(C_residual,K);
  Ug = zeros<MATTYPE>(C_residual,K);
  P = zeros<MATTYPE>(C_residual,C_residual);
  u_scale = ones<ROWVECTYPE>(N);  
  N_cb = zeros<VECTYPE>(C*B); // contains the residual space 
  N_c = zeros<VECTYPE>(C);
  Pr_b_given_c = zeros<VECTYPE>(C*B); 
  new_cell_type_rounds = zeros<VECTYPE>(C_residual+1); 
}



void dincta_full::init_cluster_cpp(unsigned Ch, float alpha) {
  cosine_normalize(Y, 0, false); // normalize columns
  
  // (2) ASSIGN CLUSTER PROBABILITIES
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  dist_mat = 2 * (1 - Y.t() * Z_cos); 
  u_scale = u * alpha + (1.- u); 
  // if C > 0, only initialize the clusters not set by the user
  // with cluster_prior
  if (Ch > 0 && Ch < K) {
    MATTYPE Rtmp = -dist_mat.rows(Ch, K-1);
    Rtmp = Rtmp.each_row()%(1./(sigma_entropy + sigma_cell_type*u_scale));
    //Rtmp *= 1./sigma;
    Rtmp.each_row() -= max(Rtmp, 0);
    Rtmp = exp(Rtmp);
    Rtmp.each_row() /= sum(Rtmp, 0);
    R.rows(Ch, K-1) = Rtmp;
  } else {
    R = -dist_mat;
    R = R.each_row()%(1./(sigma_entropy + sigma_cell_type*u_scale));
    //R *= 1/sigma;
    R.each_row() -= max(R, 0);  
    R = exp(R);
    R.each_row() /= sum(R, 0);
  }
  
  // (3) BATCH DIVERSITY STATISTICS
  
  // recompute Phi, E, O 
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      Phi.row(c*B+b) = Phi_B.row(b) % Phi_C.row(c) % u_scale;
    }
  }
  C_current_B = C_current*B;
  N_cb.head(C_current_B) = vectorise(Phi_C.head_rows(C_current).each_row() % u_scale * Phi_B.t(), 1).t(); // rows should transpose to column vector
  N_c.head(C_current) = sum(Phi_C.head_rows(C_current).each_row() % u_scale,1) + 1e-8; 
  Pr_b_given_c.head(C_current_B)  = N_cb.head(C_current_B)  / vectorise(repmat(N_c.head(C_current),1,B), 1).t(); // dim 1 concat each row
  // theta_batch_vector = theta_batch * (1 - exp(-((N_cb + 1e-6) / (nclust * tau)) ^ 2)) 
  if(tau>0){
    theta_batch_vector.head(C_current_B) = theta_batch.head(C_current_B) % ( 1.- exp(- square((N_cb.head(C_current_B) + 1e-8) / (K * tau))) );
  }
  
  E_c.head_cols(C_current) = R * (Phi_C.head_rows(C_current).t() % repmat(u_scale.t(),1,C_current)); // K x C, E K x CB
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      E.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
    }
  }
  O.head_cols(C_current_B) = R * Phi.head_rows(C_current_B).t();
  // note add the new cell type, the objective may be change, 
  // should compute one first to check the convergence. 
  compute_objective(alpha);
  objective_dincta.push_back(objective_kmeans.back());
  ran_init = true;
}

bool dincta_full::init_refine_Phi_C_and_check_u_convergence_cpp(
    float alpha, 
    bool keep_original_cell_type,
    float epsilon_cell_type_changed,
    float select_refine_fraction,
    float epsilon_changed_frequence){
  u_prior = u;
  // note that E_c for the unknown cell type is zero.
  MATTYPE T_KC = normalise(E_c.head_cols(C_current) + 1e-8, 1, 1); // Pr(c|k)
  MATTYPE Psi_C = T_KC.t() * R;
  Psi_C = normalise(Psi_C, 1, 0);
  // C - 1 not account for the unknown cell type
  //arma::urowvec cell_type_Psi_C = index_max(Psi_C, 0);
  arma::urowvec cell_type_Phi_C = index_max(Phi_C.head_rows(C_current) , 0);
  Psi_C = Phi_C.head_rows(C_current) % log( Phi_C.head_rows(C_current) / (Psi_C+1e-8));
  Psi_C.elem(find_nonfinite(Psi_C)).zeros();
  arma::rowvec kl_cell_type = sum(Psi_C, 0);
  arma::uvec kl_sort_index = sort_index(kl_cell_type, "descend");
  unsigned int top_K = ceil(select_refine_fraction * Psi_C.n_cols);
  unsigned int K_limit = as_scalar(accu(kl_cell_type>epsilon_cell_type_changed));
  if(K_limit< top_K){
    if(K_limit <10){
      top_K = 10; // avoid  k = 0
    }else{
      top_K = K_limit;
    }
  }
  unsigned int n_cell_type_c, n_cell_type_c_select, K_cell_type;
  float cell_type_select_fraction;
  arma::uvec kl_select_index = kl_sort_index.head(top_K);
  u.zeros();
  u.elem(kl_sort_index.head(top_K)).fill(1.);
  u.elem(find(kl_cell_type<=epsilon_cell_type_changed)).fill(0.);
  arma::urowvec temp;
  // not account for the unknown cell type C-1 
  // here must make sure that the last row  is the unknown cell type in Phi_C
  
  arma::uvec cell_type_ids = regspace<uvec>(0,C_known-1);
  // layout C_known, 1, C_new
  if(C_current > C_known +1){
    cell_type_ids = join_cols(cell_type_ids, regspace<uvec>(C_known+1,C_current-1));
  }
  for(unsigned int i = 0; i < cell_type_ids.size(); i++){
    unsigned int cell_type_id  = cell_type_ids(i);
    temp = cell_type_Phi_C == cell_type_id;
    n_cell_type_c = as_scalar(accu(temp));
    n_cell_type_c_select = as_scalar(accu(temp.elem(kl_select_index)));
    cell_type_select_fraction = n_cell_type_c_select*1.0 / n_cell_type_c;
    if(cell_type_select_fraction > select_refine_fraction){
      K_cell_type = floor(n_cell_type_c*select_refine_fraction);
      // find part of cell type c excluded index and set them in u to zeros in tails 
      // elem allways retrun column vector
      arma::uvec temp_select_c = temp.elem(kl_select_index.t());
      arma::uvec temp_select_c_index = kl_select_index.elem(find(temp_select_c>0));
      u.elem(temp_select_c_index.tail(n_cell_type_c_select- K_cell_type)).fill(0.);
    }
  }
  if(keep_original_cell_type){
    u.elem(find(u_original==0)).fill(0.);
  }
  
  
  if( (as_scalar(accu(abs(u_prior - u))) / (u.size() * 1.0) <=  epsilon_changed_frequence) 
        | (sum(u==1) ==0)){
    // u convergent and do not refine Phi_C
    refine_Phi_C_flag = false;
    return(true);
  } else {
    if(verbose){
      Rcout<< "refine Phi_C for " << sum(u==1) << " cells of total " << u.size() << " cells."<<endl;
    }
    compute_objective(alpha);
    objective_dincta.push_back(objective_kmeans.back());
    refine_Phi_C_flag = true;
    return(false);
  }
  
}




void dincta_full::compute_objective(float alpha) {
  float kmeans_error = as_scalar(accu(R % dist_mat)); 
  float _entropy = as_scalar(accu(safe_entropy(R) * sigma_entropy)); 
  float _cross_entropy, _kl_cell_type;
  u_scale = u * alpha + (1.- u); 
  C_current_B = C_current*B;
  _cross_entropy = as_scalar(accu((R * sigma_entropy) % ((repmat(theta_batch_vector.head(C_current_B).t(), K,1) % 
    log((O.head_cols(C_current_B) + 1) / (E.head_cols(C_current_B) + 1))) 
                                                           * Phi.head_rows(C_current_B))));
  // safe kl
  MATTYPE A = R % (
    log(R/normalise(normalise(E_c.head_cols(C_current)+1e-8,1,0)*Phi_C.head_rows(C_current), 1,0))
    % repmat(u_scale,K,1)
  );// 1 : L1 normalise, 0, for each column
  A.elem(find_nonfinite(A)).zeros(); 
  _kl_cell_type =  as_scalar(accu(A) * sigma_cell_type);
  
  objective_kmeans.push_back(kmeans_error + _entropy + _cross_entropy + _kl_cell_type );
  objective_kmeans_dist.push_back(kmeans_error);
  objective_kmeans_entropy.push_back(_entropy); 
  objective_kmeans_cross.push_back(_cross_entropy);
  objective_kmeans_kl_cell_type.push_back(_kl_cell_type);
}


bool dincta_full::check_convergence(int type) {
  float obj_new, obj_old;
  switch (type) {
  case 0: 
    // Clustering 
    // compute new window mean
    obj_old = 0;
    obj_new = 0;
    for (int i = 0; i < window_size; i++) {
      obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
      obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
    }
    if ((obj_old - obj_new) / abs(obj_old) < epsilon_kmeans) {
      return(true); 
    } else {
      return(false);
    }
  case 1:
    // dincta
    obj_old = objective_dincta[objective_dincta.size() - 2];
    obj_new = objective_dincta[objective_dincta.size() - 1];
    bool convergent_objective =(obj_old - obj_new) / abs(obj_old) < epsilon_dincta;
    bool convergent_cell_type_rounds = new_cell_type_rounds(C_new) >= min_iter_cell_type;
    if (convergent_objective & convergent_cell_type_rounds  ) {
      return(true);              
    } else {
      return(false);              
    }
  }
  
  // gives warning if we don't give default return value
  return(true);
}


int dincta_full::cluster_cpp(float alpha, int frequency_update_Phi_C) {
  int err_status = 0;
  int iter; 
  int iter_Phi_C = 0;
  Progress p(max_iter_kmeans, verbose);
  
  // Z_cos has changed
  // R has assumed to not change
  // so update Y to match new integrated data
  dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed
  for (iter = 0; iter < max_iter_kmeans; iter++) {
    p.increment();
    if (Progress::check_abort())
      return(-1);
    
    // STEP 1: Update Y
    Y = compute_Y(Z_cos, R);
    dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed
    
    // STEP 2: Update R
    err_status = update_R( alpha);
    
    if (err_status != 0) {
      // Rcout << "Compute R failed. Exiting from clustering." << endl;
      return err_status;
    }
    // STEP 3: compute objective
    compute_objective(alpha);
    if (print_loss) {
      Rcout << "epoch " << iter 
            << " kmeans loss "<< objective_kmeans.back()
            << " dist loss "<< objective_kmeans_dist.back()
            << " entropy loss "<< objective_kmeans_entropy.back()
            << " cross loss "<<  objective_kmeans_cross.back()
            << " kl cell type loss "<<  objective_kmeans_kl_cell_type.back()
            << endl;
    }
    
    
    // STEP 4: Update Phi_C
    if(update_Phi_C_flag & (iter*1.0/frequency_update_Phi_C - ceil(iter*1.0/ frequency_update_Phi_C) == 0.0) ){
      update_Phi_C(alpha);
      iter_Phi_C++;
      Phi_C_changed = true;
    }
    
    // STEP 5: check convergence
    if (iter > window_size) {
      bool convergence_status = check_convergence(0); 
      if (convergence_status) {
       if(update_Phi_C_flag& (iter*1.0/frequency_update_Phi_C - ceil(iter*1.0/ frequency_update_Phi_C) > 0.0) ){
          update_Phi_C(alpha);
          iter_Phi_C++;
          Phi_C_changed = true;
        }
        iter++;
        // Rcout << "Clustered for " << iter << " iterations" << endl;
        break;        
      }
    }
  }
  if(refine_Phi_C_flag){
    // to match the kmeans_objective value index
    kmeans_rounds.push_back(iter+1);
    refine_Phi_C_flag = false;
  }else{
    kmeans_rounds.push_back(iter);
  }
  update_Phi_C_rounds.push_back(iter_Phi_C);
  objective_dincta.push_back(objective_kmeans.back());
  new_cell_type_rounds(C_new) =  new_cell_type_rounds(C_new)+1;
  return 0;
}

int dincta_full::update_R(float alpha) { 
  //Rcout<<" begin update_R "<<endl;
  update_order = shuffle(linspace<uvec>(0, N - 1, N));
  u_scale = u * alpha + (1.- u);
  C_current_B = C_current * B;
  arma::uvec cc_list = regspace<uvec>(0, C_current - 1);
  arma::uvec ccb_list = regspace<uvec>(0, C_current_B - 1);
  
  // from the iteration logic, alpha changed only happen at each round of dincta
  // and at the end of time of each round of dincta, the Phi_C must cahnge
  if(Phi_C_changed){
    // recompute Phi, theta_batch_vector
    for(int c = 0; c < C_current; c++){
      for(unsigned b = 0; b < B; b++){
        Phi.row(c*B+b) = Phi_B.row(b) % Phi_C.row(c) % u_scale;
      }
    }
    N_cb.head(C_current_B) = vectorise(Phi_C.head_rows(C_current).each_row() % u_scale * Phi_B.t(), 1).t(); // rows should transpose to column vector
    N_c.head(C_current) = sum(Phi_C.head_rows(C_current).each_row() % u_scale,1) + 1e-8; 
    Pr_b_given_c.head(C_current_B)  = N_cb.head(C_current_B)  / vectorise(repmat(N_c.head(C_current),1,B), 1).t(); // dim 1 concat each row
    // theta_batch_vector = theta_batch * (1 - exp(-((N_cb + 1e-6) / (nclust * tau)) ^ 2)) 
    if(tau>0){
      theta_batch_vector.head(C_current_B) = theta_batch.head(C_current_B) % ( 1.- exp(- square((N_cb.head(C_current_B) + 1e-8) / (K * tau))) );
    }
    Phi_C_changed = false;
  }
  // recompute E, O
  E_c.head_cols(C_current) = R * (Phi_C.head_rows(C_current).t() % repmat(u_scale.t(),1,C_current)); // K x C, E K x CB
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      E.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
    }
  }
  O.head_cols(C_current_B) = R * Phi.head_rows(C_current_B).t();
  
  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  // randomset update 
  for (int i = 0; i < ceil(1. / block_size); i++) {    
    // gather cell updates indices
    int idx_min = i * N * block_size;
    int idx_max = min((int)((i + 1) * N * block_size), N - 1);
    // int idx_length = idx_max - idx_min + 1;
    if (idx_min > idx_max) break; // TODO: fix the loop logic so that this never happens
    uvec idx_list = linspace<uvec>(idx_min, idx_max, idx_max - idx_min + 1);
    cells_update = update_order.rows(idx_list); // can block update
    
    
    // Step 1: remove cells
    E_c.head_cols(C_current) = R.cols(cells_update) * (Phi_C(cc_list, cells_update).t() % repmat(u_scale(cells_update), 1, C_current)); // K x C, E K x CB
    for(int c = 0; c < C_current; c++){
      for(unsigned b = 0; b < B; b++){
        E_t.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
      }
    }
    
    E.head_cols(C_current_B)  -= E_t.head_cols(C_current_B) ;
    // here Phi contains scaling already
    O.head_cols(C_current_B) -= R.cols(cells_update) * Phi(ccb_list, cells_update).t();
    // this E_c is stands for N Pr(k,c)
    E_c.head_cols(C_current).zeros();
    for(int c = 0; c < C_current; c++){
      for(unsigned b = 0; b < B; b++){
        E_c.col(c) = E_c.col(c) + E.col(c*B+b);
      }
    }
    
    // Step 2: recompute R for removed cells
    R.cols(cells_update) = dist_mat.cols(cells_update) +
      +  sigma_entropy * (repmat(theta_batch_vector.head(C_current_B).t(), K,1) % log((O.head_cols(C_current_B) + 1) / (E.head_cols(C_current_B) + 1)) * Phi(ccb_list,cells_update))
      -  sigma_cell_type * 
      log( arma::max(
          normalise(normalise(E_c.head_cols(C_current)+1e-8, 1, 1)*Phi_C(cc_list, cells_update),1,0),
          1e-8 * ones(size(R.cols(cells_update)))
      )) 
      % repmat(u_scale(cells_update).t(), K, 1);
    R.cols(cells_update) = - R.cols(cells_update) /
      repmat(sigma_entropy + sigma_cell_type * u_scale.cols(cells_update), K,1);
    R.cols(cells_update) = R.cols(cells_update) - repmat(max(R.cols(cells_update), 0),K,1);
    R.cols(cells_update) = exp(R.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    R.cols(cells_update) =  (abs(R.cols(cells_update))>mu) % abs(R.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    
    // Step 3: put cells back 
    E_c.head_cols(C_current) = R.cols(cells_update) * ( Phi_C(cc_list, cells_update).t() % repmat(u_scale(cells_update), 1, C_current)) ; // K x C, E K x CB
    for(int c = 0; c < C_current; c++){
      for(unsigned b = 0; b < B; b++){
        E_t.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
      }
    }
    E.head_cols(C_current_B) += E_t.head_cols(C_current_B);
    O.head_cols(C_current_B) += R.cols(cells_update) * Phi(ccb_list, cells_update).t(); 
  }
  
  E_c.head_cols(C_current).zeros();
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      E_c.col(c) = E_c.col(c) + E.col(c*B+b);
    }
  }
  M.zeros();
  M.head_cols(C_current) = normalise(E_c.head_cols(C_current)+1e-8,1,0); // Pr(c|k) normalize coloumn
  return 0;
}

void dincta_full::update_Phi_C(float alpha) {
  //Rcout<<" begin update_Phi_C "<<endl;
  u_scale = u*alpha + (1. - u );
  VECTYPE Ru = R*u_scale.t();
  VECTYPE RI1=zeros<VECTYPE>(K) ;
  if(sum(u) < u.size() ){
    RI1 = sum(R.cols(cells_type_known_original), 1);
  }
  if(sum(RI1 < k_cluster_n_cells_outer_threshold ) >0){
    // there are new clusters in the unknown cells. 
    VECTYPE R1 = sum(R, 1);
    VECTYPE RU1 = R1 - RI1;
    arma::uvec new_clusters =  find(RI1 < k_cluster_n_cells_outer_threshold );
    
    // now use the strategy sum on the new clusters and sorted 
    // and only select neare 1 to be the main cells 
    ROWVECTYPE R_U_new_clusters_sum = sum(R(new_clusters, cells_type_unknown_original), 0);
    VECTYPE R_U_new_clusters_sum_select = R_U_new_clusters_sum(find(R_U_new_clusters_sum > new_cell_type_prob_threshold));
    int n_cells_U_new_select = sum( R_U_new_clusters_sum(find(R_U_new_clusters_sum > new_cell_type_prob_threshold)) );
    if(n_cells_U_new_select > 10){
      arma::uvec cells_type_unknown_original_select = cells_type_unknown_original(find(R_U_new_clusters_sum > new_cell_type_prob_threshold));
      // form the select cells select, based on the cells ells_type_unknown_original_select
      VECTYPE n_sample_new_clusters  = sum(R(new_clusters,cells_type_unknown_original_select),1);
      
      VECTYPE  const10 = 10* ones<VECTYPE>(n_sample_new_clusters.size());
      VECTYPE  n_main_new_clusters = max(ceil(n_sample_new_clusters*new_cell_type_main_fraction),const10);
      n_sample_new_clusters = max(floor(n_sample_new_clusters*cell_type_sample_fraction),ones<VECTYPE>(n_sample_new_clusters.size()));
      arma::uvec rare_clusters = find((n_main_new_clusters - n_sample_new_clusters)>0); 
      n_main_new_clusters(rare_clusters) =  n_sample_new_clusters(rare_clusters);
      
      // for all clusters 
      VECTYPE  n_sample_all_clusters = max(floor(R1*cell_type_sample_fraction),
                                           10* ones<VECTYPE>(R1.size()));
      n_sample_all_clusters = min(n_sample_all_clusters,ceil(R1));
      
      // gather the all cell cluster cells;
      arma::uvec sample_all_clusters_cells_gather = zeros<arma::uvec>((int)sum(n_sample_all_clusters));
      int index_sample_start=0, k=0;
      for(k = 0 ; k < n_sample_all_clusters.size(); k++){
        arma::uvec kl_sort_index = sort_index(R.row(k), "descend");
        sample_all_clusters_cells_gather.rows(index_sample_start,index_sample_start+n_sample_all_clusters(k)-1)
          = kl_sort_index.head(n_sample_all_clusters(k));
        index_sample_start += n_sample_all_clusters(k);
      }
      
      arma::uvec main_new_clusters_cells_gather = zeros<arma::uvec>((int)sum(n_main_new_clusters));
      int inc=0, index_main_start=0;
      for(inc = 0 ; inc< n_main_new_clusters.size(); inc++){
        // note that the sort index is the index in cells_type_unknown_original
        arma::uvec kl_sort_index = sort_index(R(new_clusters.at(inc)*ones<uvec>(1), cells_type_unknown_original_select), "descend");
        main_new_clusters_cells_gather.rows(index_main_start,index_main_start+n_main_new_clusters(inc)-1)
          = kl_sort_index.head(n_main_new_clusters(inc));
        index_main_start+=n_main_new_clusters(inc);
      }
      
      arma::uvec sample_all_clusters_cells = unique(sample_all_clusters_cells_gather);
      arma::uvec main_new_clusters_cells_unique = unique(main_new_clusters_cells_gather);
      
      // transform to the whole cells index
      arma::uvec main_new_clusters_cells = cells_type_unknown_original_select( main_new_clusters_cells_unique);
      
      VECTYPE R_sample_set_sum = sum( R.cols(sample_all_clusters_cells), 1);
      arma::uvec all_clusters_select = find(R_sample_set_sum>= k_cluster_n_cells_inner_threshold );
      
      
      VECTYPE R_main_set_sum = sum( R.cols(main_new_clusters_cells), 1);
      new_clusters = find(R_main_set_sum>= k_cluster_n_cells_inner_threshold );
      
      MATTYPE R_sample_all_cell_type = R(all_clusters_select, sample_all_clusters_cells);
      MATTYPE R_main_new_cell_type = R(all_clusters_select, main_new_clusters_cells);
      // compute the svd
      VECTYPE eigval;
      MATTYPE eigvec;
      VECTYPE R_sample_all_cell_type_sum  = sum(R_sample_all_cell_type,1);
      R_sample_all_cell_type_sum.elem(find(R_sample_all_cell_type_sum<1e-8)).fill(1e-8);
      MATTYPE R_sample_all_1_sqrt_inv_diag =  diagmat( 1./sqrt(R_sample_all_cell_type_sum) );
      // eigval ascending order
      MATTYPE A =  R_sample_all_1_sqrt_inv_diag * (R_sample_all_cell_type * R_sample_all_cell_type.t()) * R_sample_all_1_sqrt_inv_diag;
      eig_sym( eigval, eigvec, A);
      int C_new_current = 0, C_all_current_raw=0, C_all_current; 
      
      
      C_all_current_raw = sum( eigval > cell_type_eigval_threshold);
      if(C_all_current_raw < C_known+1){
        // if(C_all_current_raw<C_known){
        //   Rcout<<"Waring, appears "<< C_all_current_raw <<"  cell types, less than the number of known cell types " <<  C_known << endl;
        // }
        // C_all_current_raw = C_known+1;
        // here we only make sure that there must exist one cell type
        // if there is only one cell type then we will gather all the k sum to the new cell type
        if(C_all_current_raw< 1){
          C_all_current_raw = 1;
        }
      }else if(C_all_current_raw >  C_known + C_residual){
        if(verbose){
          Rcout<<"Waring, appears "<< C_all_current_raw <<" cell types, great than the all cell types " << C_known + C_residual << endl;
          Rcout<<"You can supress this waring by inrease the n.cell_type.residual to  " << C_all_current_raw - C_known << endl;
        }
        C_all_current_raw  = C_residual+C_known; // do not count the unknown cell type
      }
      
      
      
      MATTYPE U_all = find_gather_U(A,R_sample_all_cell_type_sum,C_all_current_raw,cell_type_min_cells,true,1e-8,false);
      C_all_current = U_all.n_rows;
      
      MATTYPE U_t = zeros<MATTYPE>(C_all_current, K);
      U_t.cols(all_clusters_select) = U_all;
      VECTYPE new_clusters_fraction = sum(U_t.cols(new_clusters), 1) / sum(U_all,1);
      float new_cell_type_new_clusters_fraction = 0.5;
      arma::uvec cell_types_index = sort_index(new_clusters_fraction, "descend");
      C_new_current = sum( new_clusters_fraction > new_cell_type_new_clusters_fraction );
      if(C_new_current<1){
        C_new_current = 1; 
      }
      
      MATTYPE U_new = U_all.rows(cell_types_index.head(C_new_current));
      
      
      // join cell types if centroid close too much, reduce the fake cell type
      // only use the main new cells to join cell types
      MATTYPE Phi_C_ = U_new * R_main_new_cell_type; 
      // should consider the rare cell type?
      MATTYPE Y_C_raw = Z_cos.cols(main_new_clusters_cells)*Phi_C_.t()+1e-9; // d x C_new_current
      MATTYPE Y_C,Cor_C;
      Y_C = arma::normalise(Y_C_raw,2,0);
      Cor_C = Y_C.t()* Y_C;
      Cor_C.diag().fill(-2.);
      while((Cor_C.max()>centroid_cor_threshold) & (C_new_current>=2)){
        uword im = Cor_C.index_max();
        uvec	sub = ind2sub( size(Cor_C), im);
        int c1 = sub.min(), c2 = sub.max();
        ROWVECTYPE u_combine = U_new.row(c1) + U_new.row(c2);
        VECTYPE y_raw_combine = Y_C_raw.col(c1) + Y_C_raw.col(c2);
        if(C_new_current >2){
          VECTYPE findex = ones<VECTYPE>(C_new_current);
          findex(c1) = 0;
          findex(c2) = 0;
          arma::uvec index = find(findex>0);
          U_new.head_rows(C_new_current - 2) = U_new.rows(index);
          Y_C_raw.head_cols(C_new_current - 2) = Y_C_raw.cols(index);
        }
        U_new.row(C_new_current - 2) = u_combine;
        Y_C_raw.col(C_new_current - 2) = y_raw_combine;
        C_new_current -=1;
        Y_C = arma::normalise(Y_C_raw.head_cols(C_new_current),2,0);
        Cor_C = Y_C.t()* Y_C;
        Cor_C.diag().fill(-2.);
      }
      U_new = U_new.head_rows(C_new_current);
      C_new = C_new_current;
      // layout C_known 1 C_new not need to ship the unknown column to the new positon
      C_current = C_known + C_new + 1;
      
      Phi_C.cols(main_new_clusters_cells).zeros();
     
      uvec cindex = regspace<arma::uvec>(C_known+1,C_current-1);
      Phi_C(cindex, main_new_clusters_cells) = 
        U_new * R_main_new_cell_type;
      Phi_C(cindex, main_new_clusters_cells) = normalise(Phi_C(cindex, main_new_clusters_cells), 1, 0);
      
      // // push to 1
      // arma::urowvec cell_type_position = index_max(Phi_C(cindex, main_new_clusters_cells),0);
      // ROWVECTYPE cell_type_prob = max(Phi_C(cindex, main_new_clusters_cells),0);
      // arma::uvec push_1_index = find(cell_type_prob > new_cell_type_prob_threshold);
      // cell_type_position = cell_type_position(push_1_index).t();
      // 
      // Phi_C(cindex, main_new_clusters_cells(push_1_index)).zeros();
      // for(int i =0; i<push_1_index.size();i++){
      //   Phi_C(cindex(cell_type_position(i)), main_new_clusters_cells(push_1_index(i)))=1.0;
      // }
      //Phi_C(cindex, main_new_clusters_cells) = normalise(Phi_C(cindex, main_new_clusters_cells), 1, 0);
      
      
      cells_type_known = unique(join_cols(cells_type_known_original, main_new_clusters_cells));
      VECTYPE unknown_ = ones<VECTYPE>(N);
      unknown_(cells_type_known).fill(0);
      cells_type_unknown = find(unknown_>0);
    }else{
      //Rcout <<"n_cells_U_new_select " <<n_cells_U_new_select<<"is too small. with new_cell_type_prob_threshold "<<new_cell_type_prob_threshold<<endl;
      C_new = 0;
      C_current = C_known+1;
      cells_type_known = cells_type_known_original;
      cells_type_unknown = cells_type_unknown_original;
    }
  }else{
    C_new = 0;
    C_current = C_known+1;
    cells_type_known = cells_type_known_original;
    cells_type_unknown = cells_type_unknown_original;
  }
  
  //Rcout<<" begin update_Phi_C infer "<<endl;
  Ru.elem(find(Ru< 1e-8)).fill(1e-8);
  Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) =
    (Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_known) % repmat(u_scale(cells_type_known).t(),C_current,1))
    * R.cols(cells_type_known).t() 
    % repmat(1./Ru.t(),C_current,1)
    * R.cols(cells_type_unknown);
    A_U = sqrt(alpha) * diagmat(1./sqrt(Ru)) * R.cols(cells_type_unknown);
    
    // VECTYPE eigval;
    // MATTYPE eigvec;
    // eig_sym(eigval, eigvec, A_U* A_U.t());
    // Rcout<<"eigvaule_max "<<max(eigval)<<endl;
    // float epilon_eigvalue = min(1.- max(eigval),1e-8);
    // Phi_C.cols(cells_type_unknown) = Phi_C.cols(cells_type_unknown) +
    //   Phi_C.cols(cells_type_unknown) * A_U.t() * 
    //   eigvec * arma::diagmat(1./(1. + epilon_eigvalue - eigval)) * eigvec.t() 
    //   * A_U;
    
    // VECTYPE Ru_U = R*u.t()*alpha;
    // float eigvaule_max_est = as_scalar(max(Ru_U/Ru));
    // float epilon_eigvalue = min(1.- eigvaule_max_est,1e-8);
    // Phi_C.cols(cells_type_unknown) = Phi_C.cols(cells_type_unknown) +
    //   Phi_C.cols(cells_type_unknown) * A_U.t() *  inv_sympd((1. + epilon_eigvalue) * arma::eye(K,K) - A_U* A_U.t()) * A_U;
    
    // using the solve bellow to avoid the singular problem eigvalue may near 0
    MATTYPE T_KC;
    T_KC = solve( eye(K,K) - A_U* A_U.t(), A_U*Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown).t());
    Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown)= Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) +
      T_KC.t() * A_U;
    Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) = normalise(Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown),1,0);
}


void dincta_full::moe_correct_ridge_cpp() {
  //Rcout<<" begin moe_correct_ridge_cpp"<<endl;
  Z_corr = Z_orig;
  for (int k = 0; k < K; k++) { 
    Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
    W = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
    W.row(0).zeros(); // do not remove the intercept 
    Z_corr -= W.t() * Phi_Rk;
  }
  Z_cos = arma::normalise(Z_corr, 2, 0);
}


RCPP_MODULE(dincta_full_module) {
  class_<dincta_full>("dincta_full")
  .constructor<int>()
  
  .field("Z_corr", &dincta_full::Z_corr)  
  .field("Z_orig", &dincta_full::Z_orig)  
  .field("Z_cos", &dincta_full::Z_cos)  
  .field("R", &dincta_full::R)  
  .field("Y", &dincta_full::Y)  
  .field("Phi", &dincta_full::Phi)   
  .field("Phi_B", &dincta_full::Phi_B) 
  .field("Phi_C", &dincta_full::Phi_C) 
  .field("Phi_moe", &dincta_full::Phi_moe)
  .field("Pr_b_given_c", &dincta_full::Pr_b_given_c) 
  .field("u", &dincta_full::u)
  .field("u_original", &dincta_full::u_original)
  .field("objective_kmeans", &dincta_full::objective_kmeans)
  .field("objective_kmeans_dist", &dincta_full::objective_kmeans_dist)
  .field("objective_kmeans_entropy", &dincta_full::objective_kmeans_entropy)
  .field("objective_kmeans_cross", &dincta_full::objective_kmeans_cross) 
  .field("objective_kmeans_kl_cell_type_loss", &dincta_full::objective_kmeans_kl_cell_type) 
  .field("objective_dincta", &dincta_full::objective_dincta)
  .field("centroid_max_cors", &dincta_full::centroid_max_cors)
  .field("cell_type_min_cells_vector", &dincta_full::cell_type_min_cells_vector)
  .field("dist_mat", &dincta_full::dist_mat)
  .field("ran_setup", &dincta_full::ran_setup)
  .field("ran_init", &dincta_full::ran_init)
  
  
  .field("N", &dincta_full::N)
  .field("K", &dincta_full::K)
  .field("B", &dincta_full::B)
  .field("C", &dincta_full::C)
  .field("C_known", &dincta_full::C_known)
  .field("C_new", &dincta_full::C_new)
  .field("d", &dincta_full::d)
  .field("W", &dincta_full::W)
  .field("M", &dincta_full::M)
  .field("max_iter_kmeans", &dincta_full::max_iter_kmeans)
  .field("min_iter_cell_type", &dincta_full::min_iter_cell_type) 
  //.field("frequency_update_Phi_C", &dincta_full::frequency_update_Phi_C)
  
  .field("sigma", &dincta_full::sigma)
  .field("sigma_entropy", &dincta_full::sigma_entropy)
  .field("sigma_cell_type", &dincta_full::sigma_cell_type)
  .field("theta_batch", &dincta_full::theta_batch)
  .field("lambda", &dincta_full::lambda)
  .field("mu", &dincta_full::mu)
  .field("O", &dincta_full::O) 
  .field("E", &dincta_full::E)   
  .field("E_c", &dincta_full::E_c)   
  .field("update_order", &dincta_full::update_order)    
  .field("cells_update", &dincta_full::cells_update)    
  .field("kmeans_rounds", &dincta_full::kmeans_rounds) 
  .field("new_cell_type_rounds", &dincta_full::new_cell_type_rounds) 
  .field("update_Phi_C_rounds", &dincta_full::update_Phi_C_rounds) 
  .field("epsilon_kmeans", &dincta_full::epsilon_kmeans)    
  .field("epsilon_dincta", &dincta_full::epsilon_dincta)
  
  .method("check_convergence", &dincta_full::check_convergence)
    .method("setup", &dincta_full::setup)
    .method("compute_objective", &dincta_full::compute_objective)
    .method("update_R", &dincta_full::update_R)
    .method("init_cluster_cpp", &dincta_full::init_cluster_cpp)
    .method("init_refine_Phi_C_and_check_u_convergence_cpp", 
  &dincta_full::init_refine_Phi_C_and_check_u_convergence_cpp)
  .method("cluster_cpp", &dincta_full::cluster_cpp)
  .method("moe_correct_ridge_cpp", &dincta_full::moe_correct_ridge_cpp)
  
  ;
}



// for maginal prob  distribution Pr(k,b) of Pr(k,b,c) model the cross entropy model
// useful when there are rare cell types. 

// NOTE: This is a dummy constructor, needed by Rcpp
dincta_marginal::dincta_marginal(int __K): K(__K) {}

void dincta_marginal::setup(MATTYPE& __Z, MATTYPE& __Phi, MATTYPE& __Phi_B, MATTYPE& __Phi_C,
                            MATTYPE& __Phi_moe, ROWVECTYPE __u,
                            float __sigma_entropy, 
                            float __sigma_cell_type, 
                            VECTYPE __theta_batch, // vector with size B
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
                            MATTYPE __lambda, bool __verbose, bool __print_loss) {
  
  Z_corr = MATTYPE(__Z);
  Z_orig = MATTYPE(__Z);
  Z_cos = MATTYPE(Z_orig);
  cosine_normalize(Z_cos, 0, true); // normalize columns
  
  
  Phi_B = __Phi_B;
  Phi_moe = __Phi_moe;
  u = __u;
  u_prior = 16. - u; // just make the first check_u_convergence return false
  u_original = u; 
  N = Z_corr.n_cols;
  B = Phi_B.n_rows;
  
  d = Z_corr.n_rows; 
  window_size = 3;
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_dincta = __epsilon_dincta;
  
  
  
  lambda = __lambda;
  sigma_entropy = __sigma_entropy;
  sigma_cell_type = __sigma_cell_type;
  sigma =  sigma_entropy + sigma_cell_type;
  mu = __mu; 
  sigma_prior = sigma; 
  block_size = __block_size;
  K = __K;
  C_known = __C_known;
  C_residual = __C_residual; // use for the new cell types
  k_cluster_n_cells_outer_threshold = __k_cluster_n_cells_outer_threshold;
  k_cluster_n_cells_inner_threshold = __k_cluster_n_cells_inner_threshold;
  new_cell_type_prob_threshold = __new_cell_type_prob_threshold;
  cell_type_sample_fraction = __cell_type_sample_fraction;
  new_cell_type_main_fraction = __new_cell_type_main_fraction;
  cell_type_eigval_threshold = __cell_type_eigval_threshold;
  cell_type_min_cells = __cell_type_min_cells;
  new_cell_type_min_cells = __new_cell_type_min_cells;
  centroid_cor_threshold = __centroid_cor_threshold;
  
  tau = __tau;
  // full cell type informaton 
  if(sum(u)==0){
    C_current = C_known;
    C_current_B = C_current * B;
    C_new = 0;
    C_residual = 0;
    C = C_known;
    Phi_C = __Phi_C;
    Phi = __Phi;
  }else{
    C_current = C_known + 1;
    C_current_B = C_current * B;
    C_new = 0;
    C = C_known + 1 + C_residual;
    Phi_C = zeros<MATTYPE>(C, N); 
    Phi_C.head_rows(C_current) = __Phi_C;
    Phi = zeros<MATTYPE>(C*B, N);
    Phi.head_rows(C_current*B) = __Phi;
  }
  max_iter_kmeans = __max_iter_kmeans;
  min_iter_cell_type = __min_iter_cell_type;
  verbose = __verbose;
  print_loss = __print_loss;
  
  theta_batch = __theta_batch;// unnormalize length B 
  theta_batch_vector = theta_batch; // may be normalized by N_b
  theta_batch_matrix = arma::repmat(theta_batch_vector.t(), K, 1);
  if(tau>0){
    theta_batch_vector = theta_batch % ( 1.- exp(- square((sum(Phi_B,1) + 1e-8) / (K * tau))) );
    theta_batch_matrix = arma::repmat(theta_batch_vector.t(), K, 1);
  }
  
  // unknown cell type index 
  cells_type_unknown = find(u>0); 
  cells_type_unknown_original = cells_type_unknown;
  cells_type_known = find(u==0); 
  cells_type_known_original = cells_type_known;
  update_Phi_C_flag = cells_type_unknown.size()>0;
  
  allocate_buffers();
  ran_setup = true;
  refine_Phi_C_flag = false;
}





void dincta_marginal::allocate_buffers() {
  dist_mat = zeros<MATTYPE>(K, N); 
  O = zeros<MATTYPE>(K, C*B);
  O_b = zeros<MATTYPE>(K,B);
  E = zeros<MATTYPE>(K, C*B); 
  E_b = zeros<MATTYPE>(K,B);
  E_t = zeros<MATTYPE>(K, C*B); 
  E_c = zeros<MATTYPE>(K,C);
  // todo: check the rare cell type here can get correct correction.
  W = zeros<MATTYPE>(B + 1, d); 
  Phi_Rk = zeros<MATTYPE>(B + 1, N);
  M =  zeros<MATTYPE>(K,C);
  U = zeros<MATTYPE>(C_residual,K);
  Ug = zeros<MATTYPE>(C_residual,K);
  P = zeros<MATTYPE>(C_residual,C_residual);
  u_scale = ones<ROWVECTYPE>(N);  
  N_cb = zeros<VECTYPE>(C*B); 
  N_c = zeros<VECTYPE>(C);
  Pr_b_given_c = zeros<VECTYPE>(C*B); 
  new_cell_type_rounds = zeros<VECTYPE>(C_residual+1); 
}



void dincta_marginal::init_cluster_cpp(unsigned Ch, float alpha) {
  cosine_normalize(Y, 0, false); // normalize columns
  
  // (2) ASSIGN CLUSTER PROBABILITIES
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  dist_mat = 2 * (1 - Y.t() * Z_cos); 
  u_scale = u * alpha + (1.- u);
  
  // if C > 0, only initialize the clusters not set by the user
  // with cluster_prior
  if (Ch > 0 && Ch < K) {
    MATTYPE Rtmp = -dist_mat.rows(Ch, K-1);
    // Rtmp *= 1./sigma;
    Rtmp = Rtmp.each_row()%(1./(sigma_entropy + sigma_cell_type*u_scale));
    Rtmp.each_row() -= max(Rtmp, 0);
    Rtmp = exp(Rtmp);
    Rtmp.each_row() /= sum(Rtmp, 0);
    R.rows(Ch, K-1) = Rtmp;
  } else {
    R = -dist_mat;
    R = R.each_row()%(1./(sigma_entropy + sigma_cell_type*u_scale));
    // R *= 1/sigma;
    R.each_row() -= max(R, 0);  
    R = exp(R);
    R.each_row() /= sum(R, 0);
  }
  
  // (3) BATCH DIVERSITY STATISTICS 
  // recompute Phi, E, O 
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      Phi.row(c*B+b) = Phi_B.row(b) % Phi_C.row(c) % u_scale;
    }
  }
  C_current_B = C_current*B;
  N_cb.head(C_current_B) = vectorise(Phi_C.head_rows(C_current).each_row() % u_scale * Phi_B.t(), 1).t(); // rows should transpose to column vector
  N_c.head(C_current) = sum(Phi_C.head_rows(C_current).each_row() % u_scale,1) + 1e-8; 
  Pr_b_given_c.head(C_current_B)  = N_cb.head(C_current_B)  / vectorise(repmat(N_c.head(C_current),1,B), 1).t(); // dim 1 concat each row
  
  
  E_c.head_cols(C_current) = R * (Phi_C.head_rows(C_current).t() % repmat(u_scale.t(),1,C_current)); // K x C, E K x CB
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      E.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
    }
  }
  O.head_cols(C_current_B) = R * Phi.head_rows(C_current_B).t();
  O_b.zeros();
  E_b.zeros();
  for(unsigned b = 0; b < B; b++){
    for(unsigned c = 0; c < C_current; c++){
      E_b.col(b) = E_b.col(b) + E.col(c*B+b);
      O_b.col(b) = O_b.col(b) + O.col(c*B+b);
    }
  }
  // note add the new cell type, the objective may be change, 
  // should compute one first to check the convergence. 
  compute_objective(alpha);
  objective_dincta.push_back(objective_kmeans.back());
  ran_init = true;
}

bool dincta_marginal::init_refine_Phi_C_and_check_u_convergence_cpp(
    float alpha, 
    bool keep_original_cell_type,
    float epsilon_cell_type_changed,
    float select_refine_fraction,
    float epsilon_changed_frequence){
  u_prior = u;
  // note that E_c for the unknown cell type is zero.
  MATTYPE T_KC = normalise(E_c.head_cols(C_current) + 1e-8, 1, 1); // Pr(c|k)
  MATTYPE Psi_C = T_KC.t() * R;
  Psi_C = normalise(Psi_C, 1, 0);
  // C - 1 not account for the unknown cell type
  //arma::urowvec cell_type_Psi_C = index_max(Psi_C, 0);
  arma::urowvec cell_type_Phi_C = index_max(Phi_C.head_rows(C_current) , 0);
  Psi_C = Phi_C.head_rows(C_current) % log( Phi_C.head_rows(C_current) / (Psi_C+1e-8));
  Psi_C.elem(find_nonfinite(Psi_C)).zeros();
  arma::rowvec kl_cell_type = sum(Psi_C, 0);
  arma::uvec kl_sort_index = sort_index(kl_cell_type, "descend");
  unsigned int top_K = ceil(select_refine_fraction * Psi_C.n_cols);
  unsigned int K_limit = as_scalar(accu(kl_cell_type>epsilon_cell_type_changed));
  if(K_limit< top_K){
    if(K_limit <10){
      top_K = 10; // avoid  k = 0
    }else{
      top_K = K_limit;
    }
  }
  unsigned int n_cell_type_c, n_cell_type_c_select, K_cell_type;
  float cell_type_select_fraction;
  arma::uvec kl_select_index = kl_sort_index.head(top_K);
  u.zeros();
  u.elem(kl_sort_index.head(top_K)).fill(1.);
  u.elem(find(kl_cell_type<=epsilon_cell_type_changed)).fill(0.);
  arma::urowvec temp;
  // not account for the unknown cell type C-1 
  // here must make sure that the last row  is the unknown cell type in Phi_C
  
  arma::uvec cell_type_ids = regspace<uvec>(0,C_known-1);
  // layout C_known, 1, C_new
  if(C_current > C_known +1){
    cell_type_ids = join_cols(cell_type_ids, regspace<uvec>(C_known+1,C_current-1));
  }
  for(unsigned int i = 0; i < cell_type_ids.size(); i++){
    unsigned int cell_type_id  = cell_type_ids(i);
    temp = cell_type_Phi_C == cell_type_id;
    n_cell_type_c = as_scalar(accu(temp));
    n_cell_type_c_select = as_scalar(accu(temp.elem(kl_select_index)));
    cell_type_select_fraction = n_cell_type_c_select*1.0 / n_cell_type_c;
    if(cell_type_select_fraction > select_refine_fraction){
      K_cell_type = floor(n_cell_type_c*select_refine_fraction);
      // find part of cell type c excluded index and set them in u to zeros in tails 
      // elem allways retrun column vector
      arma::uvec temp_select_c = temp.elem(kl_select_index.t());
      arma::uvec temp_select_c_index = kl_select_index.elem(find(temp_select_c>0));
      u.elem(temp_select_c_index.tail(n_cell_type_c_select- K_cell_type)).fill(0.);
    }
  }
  if(keep_original_cell_type){
    u.elem(find(u_original==0)).fill(0.);
  }
  
  
  if( (as_scalar(accu(abs(u_prior - u))) / (u.size() * 1.0) <=  epsilon_changed_frequence) 
        | (sum(u==1) ==0)){
    // u convergent and do not refine Phi_C
    refine_Phi_C_flag = false;
    return(true);
  } else {
    if(verbose){
      Rcout<< "refine Phi_C for " << sum(u==1) << " cells of total " << u.size() << " cells."<<endl;
    }
    compute_objective(alpha);
    objective_dincta.push_back(objective_kmeans.back());
    refine_Phi_C_flag = true;
    return(false);
  }
  
}


void dincta_marginal::compute_objective(float alpha) {
  float kmeans_error = as_scalar(accu(R % dist_mat)); 
  float _entropy = as_scalar(accu(safe_entropy(R) * sigma_entropy)); 
  float _cross_entropy, _kl_cell_type;
  u_scale = u * alpha + (1.- u); 
  C_current_B = C_current*B;
  _cross_entropy = as_scalar(accu((R * sigma_entropy) % ((theta_batch_matrix % log((O_b + 1) / (E_b + 1))) * Phi_B)));
  // safe kl
  MATTYPE A = R % (
    log(R/normalise(normalise(E_c.head_cols(C_current)+1e-8,1,0)*Phi_C.head_rows(C_current), 1,0))
    % repmat(u_scale,K,1)
  );// 1 : L1 normalise, 0, for each column
  A.elem(find_nonfinite(A)).zeros(); 
  _kl_cell_type =  as_scalar(accu(A) * sigma_cell_type);
  
  objective_kmeans.push_back(kmeans_error + _entropy + _cross_entropy + _kl_cell_type );
  objective_kmeans_dist.push_back(kmeans_error);
  objective_kmeans_entropy.push_back(_entropy); 
  objective_kmeans_cross.push_back(_cross_entropy);
  objective_kmeans_kl_cell_type.push_back(_kl_cell_type);
}



bool dincta_marginal::check_convergence(int type) {
  float obj_new, obj_old;
  switch (type) {
  case 0: 
    // Clustering 
    // compute new window mean
    obj_old = 0;
    obj_new = 0;
    for (int i = 0; i < window_size; i++) {
      obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
      obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
    }
    if ((obj_old - obj_new) / abs(obj_old) < epsilon_kmeans) {
      return(true); 
    } else {
      return(false);
    }
  case 1:
    // dincta
    obj_old = objective_dincta[objective_dincta.size() - 2];
    obj_new = objective_dincta[objective_dincta.size() - 1];
    bool convergent_objective =(obj_old - obj_new) / abs(obj_old) < epsilon_dincta;
    bool convergent_cell_type_rounds = new_cell_type_rounds(C_new) >= min_iter_cell_type;
    if (convergent_objective & convergent_cell_type_rounds  ) {
      return(true);              
    } else {
      return(false);              
    }
  }
  
  // gives warning if we don't give default return value
  return(true);
}



int dincta_marginal::cluster_cpp(float alpha, int frequency_update_Phi_C) {
  // Rcout<<" begin cluster_cpp "<<endl;
  int err_status = 0;
  int iter; 
  int iter_Phi_C = 0;
  Progress p(max_iter_kmeans, verbose);
  
  // Z_cos has changed
  // R has assumed to not change
  // so update Y to match new integrated data
  dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed
  for (iter = 0; iter < max_iter_kmeans; iter++) {
    p.increment();
    if (Progress::check_abort())
      return(-1);
    
    // STEP 1: Update Y
    Y = compute_Y(Z_cos, R);
    dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed
    
    // STEP 3: Update R
    err_status = update_R( alpha);
    
    if (err_status != 0) {
      // Rcout << "Compute R failed. Exiting from clustering." << endl;
      return err_status;
    }
    
    // STEP 4: compute objective
    // should in front of update Phi_C, since upadate Phi_C will change Phi_C, Phi, theta_batch 
    compute_objective(alpha);
    if (print_loss) {
      Rcout << "epoch " << iter 
            << " kmeans loss "<< objective_kmeans.back()
            << " dist loss "<< objective_kmeans_dist.back()
            << " entropy loss "<< objective_kmeans_entropy.back()
            << " cross loss "<<  objective_kmeans_cross.back()
            << " kl cell type loss "<<  objective_kmeans_kl_cell_type.back()
            << endl;
    }
    
    
    // STEP 5: Update Phi_C
    if(update_Phi_C_flag & (iter*1.0/frequency_update_Phi_C - ceil(iter*1.0/ frequency_update_Phi_C) == 0.0) ){
      update_Phi_C(alpha);
      iter_Phi_C++;
      Phi_C_changed = true;
    }
    
    // STEP 5: check convergence
    if (iter > window_size) {
      bool convergence_status = check_convergence(0); 
      if (convergence_status) {
        // update phi_C after convergent
        if(update_Phi_C_flag& (iter*1.0/frequency_update_Phi_C - ceil(iter*1.0/ frequency_update_Phi_C) > 0.0) ){
          update_Phi_C(alpha);
          iter_Phi_C++;
          Phi_C_changed = true;
        }
        iter++;
        // Rcout << "Clustered for " << iter << " iterations" << endl;
        break;        
      }
    }
  }
  if(refine_Phi_C_flag){
    // if refine, compute the objective ahead of clusterring. 
    kmeans_rounds.push_back(iter+1);
    refine_Phi_C_flag = false;
  }else{
    kmeans_rounds.push_back(iter);
  }
  update_Phi_C_rounds.push_back(iter_Phi_C);
  objective_dincta.push_back(objective_kmeans.back());
  new_cell_type_rounds(C_new) =  new_cell_type_rounds(C_new)+1;
  return 0;
}


int dincta_marginal::update_R(float alpha) { 
  // Rcout<<" begin update_R "<<endl;
  update_order = shuffle(linspace<uvec>(0, N - 1, N));
  u_scale = u * alpha + (1.- u);
  C_current_B = C_current * B;
  arma::uvec cc_list = regspace<uvec>(0, C_current - 1);
  arma::uvec ccb_list = regspace<uvec>(0, C_current_B - 1);
  
  // from the iteration logic, alpha changed only happen at each round of dincta
  // and at the end of time of each round of dincta, the Phi_C must cahnge
  if(Phi_C_changed){
    // recompute Phi
    for(int c = 0; c < C_current; c++){
      for(unsigned b = 0; b < B; b++){
        Phi.row(c*B+b) = Phi_B.row(b) % Phi_C.row(c) % u_scale;
      }
    }
    N_cb.head(C_current_B) = vectorise(Phi_C.head_rows(C_current).each_row() % u_scale * Phi_B.t(), 1).t(); // rows should transpose to column vector
    N_c.head(C_current) = sum(Phi_C.head_rows(C_current).each_row() % u_scale,1) + 1e-8; 
    Pr_b_given_c.head(C_current_B)  = N_cb.head(C_current_B)  / vectorise(repmat(N_c.head(C_current),1,B), 1).t(); // dim 1 concat each row
    Phi_C_changed = false;
  }
  // recompute  E, O 
  E_c.head_cols(C_current) = R * (Phi_C.head_rows(C_current).t() % repmat(u_scale.t(),1,C_current)); // K x C, E K x CB
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      E.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
    }
  }
  O.head_cols(C_current_B) = R * Phi.head_rows(C_current_B).t();
  
  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  // randomset update 
  for (int i = 0; i < ceil(1. / block_size); i++) {    
    // gather cell updates indices
    int idx_min = i * N * block_size;
    int idx_max = min((int)((i + 1) * N * block_size), N - 1);
    // int idx_length = idx_max - idx_min + 1;
    if (idx_min > idx_max) break; // TODO: fix the loop logic so that this never happens
    uvec idx_list = linspace<uvec>(idx_min, idx_max, idx_max - idx_min + 1);
    cells_update = update_order.rows(idx_list); // can block update
    
    
    // Step 1: remove cells
    E_c.head_cols(C_current) = R.cols(cells_update) * (Phi_C(cc_list, cells_update).t() % repmat(u_scale(cells_update), 1, C_current)); // K x C, E K x CB
    for(int c = 0; c < C_current; c++){
      for(unsigned b = 0; b < B; b++){
        E_t.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
      }
    }
    
    E.head_cols(C_current_B)  -= E_t.head_cols(C_current_B) ;
    // here Phi contains scaling already
    O.head_cols(C_current_B) -= R.cols(cells_update) * Phi(ccb_list, cells_update).t();
    
    // this E_c is stands for N Pr(k,c)
    E_c.head_cols(C_current).zeros();
    E_b.zeros();
    O_b.zeros();
    for(unsigned b = 0; b < B; b++){
      for(unsigned c = 0; c <  C_current; c++){
        E_b.col(b) = E_b.col(b) + E.col(c*B+b);
        E_c.col(c) = E_c.col(c) + E.col(c*B+b);
        O_b.col(b) = O_b.col(b) + O.col(c*B+b);
      }
    }
    
    // Step 2: recompute R for removed cells
    R.cols(cells_update) = dist_mat.cols(cells_update) +
      +  sigma_entropy * 
      (theta_batch_matrix % log((O_b + 1) / (E_b + 1)) * Phi_B.cols(cells_update))
      -  sigma_cell_type * 
        log( arma::max(
            normalise(normalise(E_c.head_cols(C_current)+1e-8, 1, 1)*Phi_C(cc_list, cells_update),1,0),
            1e-8 * ones(size(R.cols(cells_update)))
        )) 
      % repmat(u_scale(cells_update).t(), K, 1);
    R.cols(cells_update) = - R.cols(cells_update) /
      repmat(sigma_entropy + sigma_cell_type * u_scale.cols(cells_update), K,1);
    // safe exp
    R.cols(cells_update) = R.cols(cells_update) - repmat(max(R.cols(cells_update), 0),K,1);
    R.cols(cells_update) = exp(R.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    R.cols(cells_update) =  (abs(R.cols(cells_update))>mu) % abs(R.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    
    // Step 3: put cells back 
    E_c.head_cols(C_current) = R.cols(cells_update) * ( Phi_C(cc_list, cells_update).t() % repmat(u_scale(cells_update), 1, C_current)) ; // K x C, E K x CB
    for(int c = 0; c < C_current; c++){
      for(unsigned b = 0; b < B; b++){
        E_t.col(c*B+b) = E_c.col(c)*Pr_b_given_c[c*B+b];
      }
    }
    E.head_cols(C_current_B) += E_t.head_cols(C_current_B);
    O.head_cols(C_current_B) += R.cols(cells_update) * Phi(ccb_list, cells_update).t(); 
  }
  
  E_c.head_cols(C_current).zeros();
  for(int c = 0; c < C_current; c++){
    for(unsigned b = 0; b < B; b++){
      E_c.col(c) = E_c.col(c) + E.col(c*B+b);
    }
  }
  M.zeros();
  M.head_cols(C_current) = normalise(E_c.head_cols(C_current)+1e-8,1,0); // Pr(c|k) normalize coloumn
  return 0;
}


void dincta_marginal::update_Phi_C(float alpha) {
  //Rcout<<" begin update_Phi_C "<<endl;
  u_scale = u*alpha + (1. - u );
  VECTYPE Ru = R*u_scale.t();
  VECTYPE RI1=zeros<VECTYPE>(K) ;
  if(sum(u) < u.size() ){
    RI1 = sum(R.cols(cells_type_known_original), 1);
  }
  if(sum(RI1 < k_cluster_n_cells_outer_threshold ) >0){
    // there are new clusters in the unknown cells. 
    VECTYPE R1 = sum(R, 1);
    VECTYPE RU1 = R1 - RI1;
    arma::uvec new_clusters =  find(RI1 < k_cluster_n_cells_outer_threshold );
    
    // now use the strategy sum on the new clusters and sorted 
    // and only select neare 1 to be the main cells 
    ROWVECTYPE R_U_new_clusters_sum = sum(R(new_clusters, cells_type_unknown_original), 0);
    VECTYPE R_U_new_clusters_sum_select = R_U_new_clusters_sum(find(R_U_new_clusters_sum > new_cell_type_prob_threshold));
    int n_cells_U_new_select = sum( R_U_new_clusters_sum(find(R_U_new_clusters_sum > new_cell_type_prob_threshold)) );
    if(n_cells_U_new_select > 10){
      arma::uvec cells_type_unknown_original_select = cells_type_unknown_original(find(R_U_new_clusters_sum > new_cell_type_prob_threshold));
      // form the select cells select, based on the cells ells_type_unknown_original_select
      VECTYPE n_sample_new_clusters  = sum(R(new_clusters,cells_type_unknown_original_select),1);
      
      VECTYPE  const10 = 10* ones<VECTYPE>(n_sample_new_clusters.size());
      VECTYPE  n_main_new_clusters = max(ceil(n_sample_new_clusters*new_cell_type_main_fraction),const10);
      n_sample_new_clusters = max(floor(n_sample_new_clusters*cell_type_sample_fraction),ones<VECTYPE>(n_sample_new_clusters.size()));
      arma::uvec rare_clusters = find((n_main_new_clusters - n_sample_new_clusters)>0); 
      n_main_new_clusters(rare_clusters) =  n_sample_new_clusters(rare_clusters);
      
      // for all clusters 
      VECTYPE  n_sample_all_clusters = max(floor(R1*cell_type_sample_fraction),
                                           10* ones<VECTYPE>(R1.size()));
      n_sample_all_clusters = min(n_sample_all_clusters,ceil(R1));
      
      // gather the all cell cluster cells;
      arma::uvec sample_all_clusters_cells_gather = zeros<arma::uvec>((int)sum(n_sample_all_clusters));
      int index_sample_start=0, k=0;
      for(k = 0 ; k < n_sample_all_clusters.size(); k++){
        arma::uvec kl_sort_index = sort_index(R.row(k), "descend");
        sample_all_clusters_cells_gather.rows(index_sample_start,index_sample_start+n_sample_all_clusters(k)-1)
          = kl_sort_index.head(n_sample_all_clusters(k));
        index_sample_start += n_sample_all_clusters(k);
      }
      
      arma::uvec main_new_clusters_cells_gather = zeros<arma::uvec>((int)sum(n_main_new_clusters));
      int inc=0, index_main_start=0;
      for(inc = 0 ; inc< n_main_new_clusters.size(); inc++){
        // note that the sort index is the index in cells_type_unknown_original
        arma::uvec kl_sort_index = sort_index(R(new_clusters.at(inc)*ones<uvec>(1), cells_type_unknown_original_select), "descend");
        main_new_clusters_cells_gather.rows(index_main_start,index_main_start+n_main_new_clusters(inc)-1)
          = kl_sort_index.head(n_main_new_clusters(inc));
        index_main_start+=n_main_new_clusters(inc);
      }
      
      arma::uvec sample_all_clusters_cells = unique(sample_all_clusters_cells_gather);
      arma::uvec main_new_clusters_cells_unique = unique(main_new_clusters_cells_gather);
      
      // transform to the whole cells index
      arma::uvec main_new_clusters_cells = cells_type_unknown_original_select( main_new_clusters_cells_unique);
      
      VECTYPE R_sample_set_sum = sum( R.cols(sample_all_clusters_cells), 1);
      arma::uvec all_clusters_select = find(R_sample_set_sum>= k_cluster_n_cells_inner_threshold );
      
      
      VECTYPE R_main_set_sum = sum( R.cols(main_new_clusters_cells), 1);
      new_clusters = find(R_main_set_sum>= k_cluster_n_cells_inner_threshold );
      
      MATTYPE R_sample_all_cell_type = R(all_clusters_select, sample_all_clusters_cells);
      MATTYPE R_main_new_cell_type = R(all_clusters_select, main_new_clusters_cells);
      // compute the svd
      VECTYPE eigval;
      MATTYPE eigvec;
      VECTYPE R_sample_all_cell_type_sum  = sum(R_sample_all_cell_type,1);
      R_sample_all_cell_type_sum.elem(find(R_sample_all_cell_type_sum<1e-8)).fill(1e-8);
      MATTYPE R_sample_all_1_sqrt_inv_diag =  diagmat( 1./sqrt(R_sample_all_cell_type_sum) );
      // eigval ascending order
      MATTYPE A =  R_sample_all_1_sqrt_inv_diag * (R_sample_all_cell_type * R_sample_all_cell_type.t()) * R_sample_all_1_sqrt_inv_diag;
      eig_sym( eigval, eigvec, A);
      int C_new_current = 0, C_all_current_raw=0, C_all_current; 
      
      
      C_all_current_raw = sum( eigval > cell_type_eigval_threshold);
      if(C_all_current_raw < C_known+1){
        // if(C_all_current_raw<C_known){
        //   Rcout<<"Waring, appears "<< C_all_current_raw <<"  cell types, less than the number of known cell types " <<  C_known << endl;
        // }
        // C_all_current_raw = C_known+1;
        // here we only make sure that there must exist one cell type
        // if there is only one cell type then we will gather all the k sum to the new cell type
        if(C_all_current_raw< 1){
          C_all_current_raw = 1;
        }
      }else if(C_all_current_raw >  C_known + C_residual){
        if(verbose){
          Rcout<<"Waring, appears "<< C_all_current_raw <<" cell types, great than the all cell types " << C_known + C_residual << endl;
          Rcout<<"You can supress this waring by inrease the n.cell_type.residual to  " << C_all_current_raw - C_known << endl;
        }
        C_all_current_raw  = C_residual+C_known; // do not count the unknown cell type
      }
      
      
      
      MATTYPE U_all = find_gather_U(A,R_sample_all_cell_type_sum,C_all_current_raw,cell_type_min_cells,true,1e-8,false);
      C_all_current = U_all.n_rows;
      
      MATTYPE U_t = zeros<MATTYPE>(C_all_current, K);
      U_t.cols(all_clusters_select) = U_all;
      VECTYPE new_clusters_fraction = sum(U_t.cols(new_clusters), 1) / sum(U_all,1);
      float new_cell_type_new_clusters_fraction = 0.5;
      arma::uvec cell_types_index = sort_index(new_clusters_fraction, "descend");
      C_new_current = sum( new_clusters_fraction > new_cell_type_new_clusters_fraction );
      if(C_new_current<1){
        C_new_current = 1; 
      }
      
      MATTYPE U_new = U_all.rows(cell_types_index.head(C_new_current));
      
      
      // join cell types if centroid close too much, reduce the fake cell type
      // only use the main new cells to join cell types
      MATTYPE Phi_C_ = U_new * R_main_new_cell_type; 
      // should consider the rare cell type?
      MATTYPE Y_C_raw = Z_cos.cols(main_new_clusters_cells)*Phi_C_.t()+1e-9; // d x C_new_current
      MATTYPE Y_C,Cor_C;
      Y_C = arma::normalise(Y_C_raw,2,0);
      Cor_C = Y_C.t()* Y_C;
      Cor_C.diag().fill(-2.);
      while((Cor_C.max()>centroid_cor_threshold) & (C_new_current>=2)){
        uword im = Cor_C.index_max();
        uvec	sub = ind2sub( size(Cor_C), im);
        int c1 = sub.min(), c2 = sub.max();
        ROWVECTYPE u_combine = U_new.row(c1) + U_new.row(c2);
        VECTYPE y_raw_combine = Y_C_raw.col(c1) + Y_C_raw.col(c2);
        if(C_new_current >2){
          VECTYPE findex = ones<VECTYPE>(C_new_current);
          findex(c1) = 0;
          findex(c2) = 0;
          arma::uvec index = find(findex>0);
          U_new.head_rows(C_new_current - 2) = U_new.rows(index);
          Y_C_raw.head_cols(C_new_current - 2) = Y_C_raw.cols(index);
        }
        U_new.row(C_new_current - 2) = u_combine;
        Y_C_raw.col(C_new_current - 2) = y_raw_combine;
        C_new_current -=1;
        Y_C = arma::normalise(Y_C_raw.head_cols(C_new_current),2,0);
        Cor_C = Y_C.t()* Y_C;
        Cor_C.diag().fill(-2.);
      }
      U_new = U_new.head_rows(C_new_current);
      C_new = C_new_current;
      // layout C_known 1 C_new not need to ship the unknown column to the new positon
      C_current = C_known + C_new + 1;
      
      Phi_C.cols(main_new_clusters_cells).zeros();
      
      uvec cindex = regspace<arma::uvec>(C_known+1,C_current-1);
      Phi_C(cindex, main_new_clusters_cells) = 
        U_new * R_main_new_cell_type;
      Phi_C(cindex, main_new_clusters_cells) = normalise(Phi_C(cindex, main_new_clusters_cells), 1, 0);
      
      // // push to 1
      // arma::urowvec cell_type_position = index_max(Phi_C(cindex, main_new_clusters_cells),0);
      // ROWVECTYPE cell_type_prob = max(Phi_C(cindex, main_new_clusters_cells),0);
      // arma::uvec push_1_index = find(cell_type_prob > new_cell_type_prob_threshold);
      // cell_type_position = cell_type_position(push_1_index).t();
      // 
      // Phi_C(cindex, main_new_clusters_cells(push_1_index)).zeros();
      // for(int i =0; i<push_1_index.size();i++){
      //   Phi_C(cindex(cell_type_position(i)), main_new_clusters_cells(push_1_index(i)))=1.0;
      // }
      //Phi_C(cindex, main_new_clusters_cells) = normalise(Phi_C(cindex, main_new_clusters_cells), 1, 0);
      
      
      cells_type_known = unique(join_cols(cells_type_known_original, main_new_clusters_cells));
      VECTYPE unknown_ = ones<VECTYPE>(N);
      unknown_(cells_type_known).fill(0);
      cells_type_unknown = find(unknown_>0);
    }else{
      //Rcout <<"n_cells_U_new_select " <<n_cells_U_new_select<<"is too small. with new_cell_type_prob_threshold "<<new_cell_type_prob_threshold<<endl;
      C_new = 0;
      C_current = C_known+1;
      cells_type_known = cells_type_known_original;
      cells_type_unknown = cells_type_unknown_original;
    }
  }else{
    C_new = 0;
    C_current = C_known+1;
    cells_type_known = cells_type_known_original;
    cells_type_unknown = cells_type_unknown_original;
  }
  
  //Rcout<<" begin update_Phi_C infer "<<endl;
  Ru.elem(find(Ru< 1e-8)).fill(1e-8);
  Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) =
    (Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_known) % repmat(u_scale(cells_type_known).t(),C_current,1))
    * R.cols(cells_type_known).t() 
    % repmat(1./Ru.t(),C_current,1)
    * R.cols(cells_type_unknown);
    A_U = sqrt(alpha) * diagmat(1./sqrt(Ru)) * R.cols(cells_type_unknown);
    
    // VECTYPE eigval;
    // MATTYPE eigvec;
    // eig_sym(eigval, eigvec, A_U* A_U.t());
    // Rcout<<"eigvaule_max "<<max(eigval)<<endl;
    // float epilon_eigvalue = min(1.- max(eigval),1e-8);
    // Phi_C.cols(cells_type_unknown) = Phi_C.cols(cells_type_unknown) +
    //   Phi_C.cols(cells_type_unknown) * A_U.t() * 
    //   eigvec * arma::diagmat(1./(1. + epilon_eigvalue - eigval)) * eigvec.t() 
    //   * A_U;
    
    // VECTYPE Ru_U = R*u.t()*alpha;
    // float eigvaule_max_est = as_scalar(max(Ru_U/Ru));
    // float epilon_eigvalue = min(1.- eigvaule_max_est,1e-8);
    // Phi_C.cols(cells_type_unknown) = Phi_C.cols(cells_type_unknown) +
    //   Phi_C.cols(cells_type_unknown) * A_U.t() *  inv_sympd((1. + epilon_eigvalue) * arma::eye(K,K) - A_U* A_U.t()) * A_U;
    
    // using the solve bellow to avoid the singular problem eigvalue may near 0
    MATTYPE T_KC;
    T_KC = solve( eye(K,K) - A_U* A_U.t(), A_U*Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown).t());
    Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown)= Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) +
      T_KC.t() * A_U;
    Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) = normalise(Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown),1,0);
}


void dincta_marginal::moe_correct_ridge_cpp() {
  //Rcout<<" begin moe_correct_ridge_cpp"<<endl;
  Z_corr = Z_orig;
  for (int k = 0; k < K; k++) { 
    Phi_Rk = Phi_moe * arma::diagmat(R.row(k));
    W = arma::inv(Phi_Rk * Phi_moe.t() + lambda) * Phi_Rk * Z_orig.t();
    W.row(0).zeros(); // do not remove the intercept 
    Z_corr -= W.t() * Phi_Rk;
  }
  Z_cos = arma::normalise(Z_corr, 2, 0);
}


RCPP_MODULE(dincta_marginal_module) {
  class_<dincta_marginal>("dincta_marginal")
  .constructor<int>()
  
  .field("Z_corr", &dincta_marginal::Z_corr)
  .field("Z_orig", &dincta_marginal::Z_orig)  
  .field("Z_cos", &dincta_marginal::Z_cos)  
  .field("R", &dincta_marginal::R)  
  .field("Y", &dincta_marginal::Y)  
  .field("Phi", &dincta_marginal::Phi)   
  .field("Phi_B", &dincta_marginal::Phi_B) 
  .field("Phi_C", &dincta_marginal::Phi_C) 
  .field("Phi_moe", &dincta_marginal::Phi_moe)
  .field("Pr_b_given_c", &dincta_marginal::Pr_b_given_c) 
  .field("u", &dincta_marginal::u)
  .field("objective_kmeans", &dincta_marginal::objective_kmeans)
  .field("objective_kmeans_dist", &dincta_marginal::objective_kmeans_dist)
  .field("objective_kmeans_entropy", &dincta_marginal::objective_kmeans_entropy)
  .field("objective_kmeans_cross", &dincta_marginal::objective_kmeans_cross) 
  .field("objective_kmeans_kl_cell_type_loss", &dincta_marginal::objective_kmeans_kl_cell_type) 
  .field("objective_dincta", &dincta_marginal::objective_dincta)
  .field("centroid_max_cors", &dincta_marginal::centroid_max_cors)
  .field("cell_type_min_cells_vector", &dincta_marginal::cell_type_min_cells_vector)
  .field("dist_mat", &dincta_marginal::dist_mat)
  .field("ran_setup", &dincta_marginal::ran_setup)
  .field("ran_init", &dincta_marginal::ran_init)
  
  
  .field("N", &dincta_marginal::N)
  .field("K", &dincta_marginal::K)
  .field("B", &dincta_marginal::B)
  .field("C", &dincta_marginal::C)
  .field("C_known", &dincta_marginal::C_known)
  .field("C_new", &dincta_marginal::C_new)
  .field("d", &dincta_marginal::d)
  .field("W", &dincta_marginal::W)
  .field("M", &dincta_marginal::M)
  .field("max_iter_kmeans", &dincta_marginal::max_iter_kmeans)
  //.field("frequency_update_Phi_C", &dincta_marginal::frequency_update_Phi_C)
  
  .field("sigma", &dincta_marginal::sigma)
  .field("sigma_entropy", &dincta_marginal::sigma_entropy)
  .field("sigma_cell_type", &dincta_marginal::sigma_cell_type)
  .field("theta_batch", &dincta_marginal::theta_batch)
  .field("lambda", &dincta_marginal::lambda)
  .field("mu", &dincta_marginal::mu)
  .field("u_original", &dincta_marginal::u_original)
  .field("O", &dincta_marginal::O) 
  .field("E", &dincta_marginal::E)   
  .field("E_c", &dincta_marginal::E_c)   
  .field("update_order", &dincta_marginal::update_order)    
  .field("cells_update", &dincta_marginal::cells_update)    
  .field("kmeans_rounds", &dincta_marginal::kmeans_rounds)    
  .field("update_Phi_C_rounds", &dincta_marginal::update_Phi_C_rounds) 
  .field("new_cell_type_rounds", &dincta_marginal::new_cell_type_rounds) 
  .field("epsilon_kmeans", &dincta_marginal::epsilon_kmeans)    
  .field("epsilon_dincta", &dincta_marginal::epsilon_dincta)
  
  //.method("init_cluster", &dincta_marginal::init_cluster_cpp)
    .method("check_convergence", &dincta_marginal::check_convergence)
    .method("setup", &dincta_marginal::setup)
    .method("compute_objective", &dincta_marginal::compute_objective)
    .method("update_R", &dincta_marginal::update_R)
    .method("init_cluster_cpp", &dincta_marginal::init_cluster_cpp)
    .method("init_refine_Phi_C_and_check_u_convergence_cpp", 
  &dincta_marginal::init_refine_Phi_C_and_check_u_convergence_cpp)
  .method("cluster_cpp", &dincta_marginal::cluster_cpp)
  .method("moe_correct_ridge_cpp", &dincta_marginal::moe_correct_ridge_cpp)
  
  ;
}






// Soft kmeans classification  model  

// NOTE: This is a dummy constructor, needed by Rcpp
dincta_classification::dincta_classification(int __K): K(__K) {}



void dincta_classification::setup(MATTYPE& __Z, MATTYPE& __Phi_C,
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
                                  bool __verbose, bool __print_loss) {
  
  Z_corr = MATTYPE(__Z);
  Z_orig = MATTYPE(__Z);
  Z_cos = MATTYPE(Z_orig);
  cosine_normalize(Z_cos, 0, true); // normalize columns
  
  u = __u;
  u_prior = 16. - u; // just make the first check_u_convergence return false
  u_original = u; 
  N = Z_corr.n_cols;
  
  d = Z_corr.n_rows; 
  window_size = 3;
  epsilon_kmeans = __epsilon_kmeans;
  epsilon_dincta = __epsilon_dincta;
  
  
  
  sigma_entropy = __sigma_entropy;
  sigma_cell_type = __sigma_cell_type;
  sigma =  sigma_entropy + sigma_cell_type;
  mu = __mu; 
  sigma_prior = sigma; // where use ? 
  block_size = __block_size;
  K = __K;
  C_known = __C_known;
  C_residual = __C_residual; // use for the new cell types
  
  k_cluster_n_cells_outer_threshold = __k_cluster_n_cells_outer_threshold;
  k_cluster_n_cells_inner_threshold = __k_cluster_n_cells_inner_threshold;
  new_cell_type_prob_threshold = __new_cell_type_prob_threshold;
  cell_type_sample_fraction = __cell_type_sample_fraction;
  new_cell_type_main_fraction = __new_cell_type_main_fraction;
  cell_type_eigval_threshold = __cell_type_eigval_threshold;
  cell_type_min_cells = __cell_type_min_cells;
  new_cell_type_min_cells = __new_cell_type_min_cells;
  centroid_cor_threshold = __centroid_cor_threshold;
  // full cell type informaton 
  if(sum(u)==0){
    C_current = C_known;
    C_new = 0;
    C_residual = 0;
    C = C_known;
    Phi_C = __Phi_C;
  }else{
    C_current = C_known + 1;
    C_new = 0;
    C = C_known + 1 + C_residual;
    Phi_C = zeros<MATTYPE>(C, N); 
    Phi_C.head_rows(C_current) = __Phi_C;
  }
  max_iter_kmeans = __max_iter_kmeans;
  // min_iter_cell_type = __min_iter_cell_type;
  verbose = __verbose;
  print_loss = __print_loss;
  
  // unknown cell type index 
  cells_type_unknown = find(u>0); 
  cells_type_unknown_original = cells_type_unknown;
  cells_type_known = find(u==0); 
  cells_type_known_original = cells_type_known;
  update_Phi_C_flag = cells_type_unknown.size()>0;
  
  allocate_buffers();
  ran_setup = true;
  refine_Phi_C_flag = false;
}





void dincta_classification::allocate_buffers() {
  dist_mat = zeros<MATTYPE>(K, N); 
  O = zeros<MATTYPE>(K, C);
  M =  zeros<MATTYPE>(K,C);
  U = zeros<MATTYPE>(C_residual,K);
  Ug = zeros<MATTYPE>(C_residual,K);
  P = zeros<MATTYPE>(C_residual,C_residual);
  u_scale = ones<ROWVECTYPE>(N);  
  N_c = zeros<VECTYPE>(C);
}



void dincta_classification::init_cluster_cpp(unsigned Ch, float alpha) {
  //Rcout<<" begin init_cluster_cpp "<<endl;
  cosine_normalize(Y, 0, false); // normalize columns
  
  // (2) ASSIGN CLUSTER PROBABILITIES
  // using a nice property of cosine distance,
  // compute squared distance directly with cross product
  dist_mat = 2 * (1 - Y.t() * Z_cos); 
  u_scale = u * alpha + (1.- u); 
  
  // if C > 0, only initialize the clusters not set by the user
  // with cluster_prior
  if (Ch > 0 && Ch < K) {
    MATTYPE Rtmp = -dist_mat.rows(Ch, K-1);
    Rtmp = Rtmp.each_row()%(1./(sigma_entropy + sigma_cell_type*u_scale));
    //Rtmp *= 1./sigma;
    Rtmp.each_row() -= max(Rtmp, 0);
    Rtmp = exp(Rtmp);
    Rtmp.each_row() /= sum(Rtmp, 0);
    R.rows(Ch, K-1) = Rtmp;
  } else {
    R = -dist_mat;
    R = R.each_row()%(1./(sigma_entropy + sigma_cell_type*u_scale));
    //R *= 1/sigma;
    R.each_row() -= max(R, 0);  
    R = exp(R);
    R.each_row() /= sum(R, 0);
  }
  
  // (3) BATCH DIVERSITY STATISTICS
  // recompute Phi, E, O 
  N_c.head(C_current) = sum(Phi_C.head_rows(C_current).each_row() % u_scale,1) + 1e-8; 
  O.head_cols(C_current) = R * (Phi_C.head_rows(C_current).t() % repmat(u_scale.t(),1,C_current)); // K x C, E K x CB
  
  // note add the new cell type, the objective may be change, 
  // should compute one first to check the convergence. 
  compute_objective(alpha);
  objective_dincta.push_back(objective_kmeans.back());
  ran_init = true;
}

bool dincta_classification::init_refine_Phi_C_and_check_u_convergence_cpp(
    float alpha, 
    bool keep_original_cell_type,
    float epsilon_cell_type_changed,
    float select_refine_fraction,
    float epsilon_changed_frequence){
  //Rcout<<" begin init_refine_Phi_C_and_check_u_convergence_cpp "<<endl;
  u_prior = u;
  // note that E_c for the unknown cell type is zero.
  MATTYPE T_KC = normalise(O.head_cols(C_current) + 1e-8, 1, 1); // Pr(c|k)
  MATTYPE Psi_C = T_KC.t() * R;
  Psi_C = normalise(Psi_C, 1, 0);
  arma::urowvec cell_type_Phi_C = index_max(Phi_C.head_rows(C_current) , 0);
  Psi_C = Phi_C.head_rows(C_current) % log( Phi_C.head_rows(C_current) / (Psi_C+1e-8));
  Psi_C.elem(find_nonfinite(Psi_C)).zeros();
  arma::rowvec kl_cell_type = sum(Psi_C, 0);
  arma::uvec kl_sort_index = sort_index(kl_cell_type, "descend");
  unsigned int top_K = ceil(select_refine_fraction * Psi_C.n_cols);
  unsigned int K_limit = as_scalar(accu(kl_cell_type>epsilon_cell_type_changed));
  if(K_limit< top_K){
    if(K_limit <10){
      top_K = 10; // avoid  k = 0
    }else{
      top_K = K_limit;
    }
  }
  unsigned int n_cell_type_c, n_cell_type_c_select, K_cell_type;
  float cell_type_select_fraction;
  arma::uvec kl_select_index = kl_sort_index.head(top_K);
  u.zeros();
  u.elem(kl_sort_index.head(top_K)).fill(1.);
  u.elem(find(kl_cell_type<=epsilon_cell_type_changed)).fill(0.);
  arma::urowvec temp;
  // not account for the unknown cell type C-1 
  // here must make sure that the last row  is the unknown cell type in Phi_C
  
  arma::uvec cell_type_ids = regspace<uvec>(0,C_known-1);
  // layout C_known, 1, C_new
  if(C_current > C_known +1){
    cell_type_ids = join_cols(cell_type_ids, regspace<uvec>(C_known+1,C_current-1));
  }
  for(unsigned int i = 0; i < cell_type_ids.size(); i++){
    unsigned int cell_type_id  = cell_type_ids(i);
    temp = cell_type_Phi_C == cell_type_id;
    n_cell_type_c = as_scalar(accu(temp));
    n_cell_type_c_select = as_scalar(accu(temp.elem(kl_select_index)));
    cell_type_select_fraction = n_cell_type_c_select*1.0 / n_cell_type_c;
    if(cell_type_select_fraction > select_refine_fraction){
      K_cell_type = floor(n_cell_type_c*select_refine_fraction);
      // find part of cell type c excluded index and set them in u to zeros in tails 
      // elem allways retrun column vector
      arma::uvec temp_select_c = temp.elem(kl_select_index.t());
      arma::uvec temp_select_c_index = kl_select_index.elem(find(temp_select_c>0));
      u.elem(temp_select_c_index.tail(n_cell_type_c_select- K_cell_type)).fill(0.);
    }
  }
  if(keep_original_cell_type){
    u.elem(find(u_original==0)).fill(0.);
  }
  
  
  // compare u and u_prior to decide refine or not
  if( (as_scalar(accu(abs(u_prior - u))) / (u.size() * 1.0) <=  epsilon_changed_frequence) 
        | (sum(u==1) ==0)){
    // u convergent and do not refine Phi_C
    refine_Phi_C_flag = false;
    return(true);
  } else {
    if(verbose){
      Rcout<< "refine Phi_C for " << sum(u==1) << " cells of total " << u.size() << " cells."<<endl;
    }
    compute_objective(alpha);
    objective_dincta.push_back(objective_kmeans.back());
    refine_Phi_C_flag = true;
    return(false);
  }
  
}


void dincta_classification::compute_objective(float alpha) {
  float kmeans_error = as_scalar(accu(R % dist_mat)); 
  float _entropy = as_scalar(accu(safe_entropy(R) * sigma_entropy)); 
  float _kl_cell_type;
  u_scale = u * alpha + (1.- u); 
  // safe kl
  MATTYPE A = R % (
    log(R/normalise(normalise(O.head_cols(C_current)+1e-8,1,0)*Phi_C.head_rows(C_current), 1,0))
    % repmat(u_scale,K,1)
  );// 1 : L1 normalise, 0, for each column
  A.elem(find_nonfinite(A)).zeros(); 
  _kl_cell_type =  as_scalar(accu(A) * sigma_cell_type);
  
  objective_kmeans.push_back(kmeans_error + _entropy + _kl_cell_type );
  objective_kmeans_dist.push_back(kmeans_error);
  objective_kmeans_entropy.push_back(_entropy); 
  objective_kmeans_kl_cell_type.push_back(_kl_cell_type);
}



bool dincta_classification::check_convergence(int type) {
  float obj_new, obj_old;
  switch (type) {
  case 0: 
    // Clustering 
    // compute new window mean
    obj_old = 0;
    obj_new = 0;
    for (int i = 0; i < window_size; i++) {
      obj_old += objective_kmeans[objective_kmeans.size() - 2 - i];
      obj_new += objective_kmeans[objective_kmeans.size() - 1 - i];
    }
    if ((obj_old - obj_new) / abs(obj_old) < epsilon_kmeans) {
      return(true); 
    } else {
      return(false);
    }
  case 1:
    // dincta used for refine?
    obj_old = objective_dincta[objective_dincta.size() - 2];
    obj_new = objective_dincta[objective_dincta.size() - 1];
    if ((obj_old - obj_new) / abs(obj_old) < epsilon_dincta) {
      return(true);              
    } else {
      return(false);              
    }
  }
  
  // gives warning if we don't give default return value
  return(true);
}



int dincta_classification::cluster_cpp(float alpha, int frequency_update_Phi_C) {
  //Rcout<<" begin cluster_cpp "<<endl;
  int err_status = 0;
  int iter; 
  int iter_Phi_C = 0;
  Progress p(max_iter_kmeans, verbose);
  
  // Z_cos has changed
  // R has assumed to not change
  // so update Y to match new integrated data
  dist_mat = 2 * (1 - Y.t() * Z_cos); // Z_cos was changed
  for (iter = 0; iter < max_iter_kmeans; iter++) {
    p.increment();
    if (Progress::check_abort())
      return(-1);
    
    // STEP 1: Update Y
    Y = compute_Y(Z_cos, R);
    dist_mat = 2 * (1 - Y.t() * Z_cos); // Y was changed
    
    // STEP 3: Update R
    err_status = update_R( alpha);
    
    if (err_status != 0) {
      // Rcout << "Compute R failed. Exiting from clustering." << endl;
      return err_status;
    }
    
    // STEP 4: compute objective
    // should in front of update Phi_C, since upadate Phi_C will change Phi_C, Phi, theta_batch 
    compute_objective(alpha);
    if (print_loss) {
      Rcout << "epoch " << iter 
            << " kmeans loss "<< objective_kmeans.back()
            << " dist loss "<< objective_kmeans_dist.back()
            << " entropy loss "<< objective_kmeans_entropy.back()
            << " kl cell type loss "<<  objective_kmeans_kl_cell_type.back()
            << endl;
    }
    
    
    // STEP 5: Update Phi_C
    if(update_Phi_C_flag & (iter*1.0/frequency_update_Phi_C - ceil(iter*1.0/ frequency_update_Phi_C) == 0.0) ){
      update_Phi_C(alpha);
      iter_Phi_C++;
      Phi_C_changed = true;
    }
    
    // STEP 5: check convergence
    if (iter > window_size) {
      bool convergence_status = check_convergence(0); 
      if (convergence_status) {
        if(update_Phi_C_flag& (iter*1.0/frequency_update_Phi_C - ceil(iter*1.0/ frequency_update_Phi_C) > 0.0) ){
          update_Phi_C(alpha);
          iter_Phi_C++;
          Phi_C_changed = true;
        }
        iter++;
        // Rcout << "Clustered for " << iter << " iterations" << endl;
        break;        
      }
    }
  }
  M.zeros();
  M.head_cols(C_current) = normalise(O.head_cols(C_current)+1e-8,1,0); // Pr(c|k) normalize coloumn
  if(refine_Phi_C_flag){
    // to match the kmeans_objective value index
    kmeans_rounds.push_back(iter+1);
    refine_Phi_C_flag = false;
  }else{
    kmeans_rounds.push_back(iter);
  }
  update_Phi_C_rounds.push_back(iter_Phi_C);
  objective_dincta.push_back(objective_kmeans.back());
  return 0;
}


int dincta_classification::update_R(float alpha) { 
  //Rcout<<" begin update_R "<<endl;
  update_order = shuffle(linspace<uvec>(0, N - 1, N));
  u_scale = u * alpha + (1.- u);
  arma::uvec cc_list = regspace<uvec>(0, C_current - 1);
  
  // from the iteration logic, alpha changed only happen at each round of dincta
  // and at the end of time of each round of dincta, the Phi_C must cahnge
  // recompute Phi, E, O 
  if(Phi_C_changed){
    N_c.head(C_current) = sum(Phi_C.head_rows(C_current).each_row() % u_scale,1) + 1e-8; 
    Phi_C_changed = false;
  }
  O.head_cols(C_current) = R * (Phi_C.head_rows(C_current).t() % repmat(u_scale.t(),1,C_current)); // K x C, E K x CB
  // GENERAL CASE: online updates, in blocks of size (N * block_size)
  // randomset update 
  for (int i = 0; i < ceil(1. / block_size); i++) {    
    // gather cell updates indices
    int idx_min = i * N * block_size;
    int idx_max = min((int)((i + 1) * N * block_size), N - 1);
    // int idx_length = idx_max - idx_min + 1;
    if (idx_min > idx_max) break; // TODO: fix the loop logic so that this never happens
    uvec idx_list = linspace<uvec>(idx_min, idx_max, idx_max - idx_min + 1);
    cells_update = update_order.rows(idx_list); // can block update
    
    
    // Step 1: remove cells
    O.head_cols(C_current) -=  R.cols(cells_update) * (Phi_C(cc_list, cells_update).t() % repmat(u_scale(cells_update), 1, C_current)); 
    
    // Step 2: recompute R for removed cells
    R.cols(cells_update) = dist_mat.cols(cells_update) +
      -  sigma_cell_type * 
      log( arma::max(
          normalise(normalise(O.head_cols(C_current)+1e-8, 1, 1)*Phi_C(cc_list, cells_update),1,0),
          1e-8 * ones(size(R.cols(cells_update)))
      )) 
      % repmat(u_scale(cells_update).t(), K, 1);
    R.cols(cells_update) = - R.cols(cells_update) /
      repmat(sigma_entropy + sigma_cell_type * u_scale.cols(cells_update), K,1);
    R.cols(cells_update) = R.cols(cells_update) - repmat(max(R.cols(cells_update), 0),K,1);
    R.cols(cells_update) = exp(R.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    R.cols(cells_update) =  (abs(R.cols(cells_update))>mu) % abs(R.cols(cells_update));
    R.cols(cells_update) = normalise(R.cols(cells_update), 1, 0); // L1 norm columns
    
    // Step 3: put cells back 
    O.head_cols(C_current) += R.cols(cells_update) * ( Phi_C(cc_list, cells_update).t() % repmat(u_scale(cells_update), 1, C_current)) ; // K x C, E K x CB
  }
  return 0;
}


void dincta_classification::update_Phi_C(float alpha) {
  //Rcout<<" begin update_Phi_C "<<endl;
  u_scale = u*alpha + (1. - u );
  VECTYPE Ru = R*u_scale.t();
  VECTYPE RI1=zeros<VECTYPE>(K) ;
  if(sum(u) < u.size() ){
    RI1 = sum(R.cols(cells_type_known_original), 1);
  }
  if(sum(RI1 < k_cluster_n_cells_outer_threshold ) >0){
    // there are new clusters in the unknown cells. 
    VECTYPE R1 = sum(R, 1);
    VECTYPE RU1 = R1 - RI1;
    arma::uvec new_clusters =  find(RI1 < k_cluster_n_cells_outer_threshold );
    
    // now use the strategy sum on the new clusters and sorted 
    // and only select neare 1 to be the main cells 
    ROWVECTYPE R_U_new_clusters_sum = sum(R(new_clusters, cells_type_unknown_original), 0);
    VECTYPE R_U_new_clusters_sum_select = R_U_new_clusters_sum(find(R_U_new_clusters_sum > new_cell_type_prob_threshold));
    int n_cells_U_new_select = sum( R_U_new_clusters_sum(find(R_U_new_clusters_sum > new_cell_type_prob_threshold)) );
    if(n_cells_U_new_select > 10){
      arma::uvec cells_type_unknown_original_select = cells_type_unknown_original(find(R_U_new_clusters_sum > new_cell_type_prob_threshold));
      // form the select cells select, based on the cells ells_type_unknown_original_select
      VECTYPE n_sample_new_clusters  = sum(R(new_clusters,cells_type_unknown_original_select),1);
      
      VECTYPE  const10 = 10* ones<VECTYPE>(n_sample_new_clusters.size());
      VECTYPE  n_main_new_clusters = max(ceil(n_sample_new_clusters*new_cell_type_main_fraction),const10);
      n_sample_new_clusters = max(floor(n_sample_new_clusters*cell_type_sample_fraction),ones<VECTYPE>(n_sample_new_clusters.size()));
      arma::uvec rare_clusters = find((n_main_new_clusters - n_sample_new_clusters)>0); 
      n_main_new_clusters(rare_clusters) =  n_sample_new_clusters(rare_clusters);
      
      // for all clusters 
      VECTYPE  n_sample_all_clusters = max(floor(R1*cell_type_sample_fraction),
                                           10* ones<VECTYPE>(R1.size()));
      n_sample_all_clusters = min(n_sample_all_clusters,ceil(R1));
      
      // gather the all cell cluster cells;
      arma::uvec sample_all_clusters_cells_gather = zeros<arma::uvec>((int)sum(n_sample_all_clusters));
      int index_sample_start=0, k=0;
      for(k = 0 ; k < n_sample_all_clusters.size(); k++){
        arma::uvec kl_sort_index = sort_index(R.row(k), "descend");
        sample_all_clusters_cells_gather.rows(index_sample_start,index_sample_start+n_sample_all_clusters(k)-1)
          = kl_sort_index.head(n_sample_all_clusters(k));
        index_sample_start += n_sample_all_clusters(k);
      }
      
      arma::uvec main_new_clusters_cells_gather = zeros<arma::uvec>((int)sum(n_main_new_clusters));
      int inc=0, index_main_start=0;
      for(inc = 0 ; inc< n_main_new_clusters.size(); inc++){
        // note that the sort index is the index in cells_type_unknown_original
        arma::uvec kl_sort_index = sort_index(R(new_clusters.at(inc)*ones<uvec>(1), cells_type_unknown_original_select), "descend");
        main_new_clusters_cells_gather.rows(index_main_start,index_main_start+n_main_new_clusters(inc)-1)
          = kl_sort_index.head(n_main_new_clusters(inc));
        index_main_start+=n_main_new_clusters(inc);
      }
      
      arma::uvec sample_all_clusters_cells = unique(sample_all_clusters_cells_gather);
      arma::uvec main_new_clusters_cells_unique = unique(main_new_clusters_cells_gather);
      
      // transform to the whole cells index
      arma::uvec main_new_clusters_cells = cells_type_unknown_original_select( main_new_clusters_cells_unique);
      
      VECTYPE R_sample_set_sum = sum( R.cols(sample_all_clusters_cells), 1);
      arma::uvec all_clusters_select = find(R_sample_set_sum>= k_cluster_n_cells_inner_threshold );
      
      
      VECTYPE R_main_set_sum = sum( R.cols(main_new_clusters_cells), 1);
      new_clusters = find(R_main_set_sum>= k_cluster_n_cells_inner_threshold );
      
      MATTYPE R_sample_all_cell_type = R(all_clusters_select, sample_all_clusters_cells);
      MATTYPE R_main_new_cell_type = R(all_clusters_select, main_new_clusters_cells);
      // compute the svd
      VECTYPE eigval;
      MATTYPE eigvec;
      VECTYPE R_sample_all_cell_type_sum  = sum(R_sample_all_cell_type,1);
      R_sample_all_cell_type_sum.elem(find(R_sample_all_cell_type_sum<1e-8)).fill(1e-8);
      MATTYPE R_sample_all_1_sqrt_inv_diag =  diagmat( 1./sqrt(R_sample_all_cell_type_sum) );
      // eigval ascending order
      MATTYPE A =  R_sample_all_1_sqrt_inv_diag * (R_sample_all_cell_type * R_sample_all_cell_type.t()) * R_sample_all_1_sqrt_inv_diag;
      eig_sym( eigval, eigvec, A);
      int C_new_current = 0, C_all_current_raw=0, C_all_current; 
      
      
      C_all_current_raw = sum( eigval > cell_type_eigval_threshold);
      if(C_all_current_raw < C_known+1){
        // if(C_all_current_raw<C_known){
        //   Rcout<<"Waring, appears "<< C_all_current_raw <<"  cell types, less than the number of known cell types " <<  C_known << endl;
        // }
        // C_all_current_raw = C_known+1;
        // here we only make sure that there must exist one cell type
        // if there is only one cell type then we will gather all the k sum to the new cell type
        if(C_all_current_raw< 1){
          C_all_current_raw = 1;
        }
      }else if(C_all_current_raw >  C_known + C_residual){
        if(verbose){
          Rcout<<"Waring, appears "<< C_all_current_raw <<" cell types, great than the all cell types " << C_known + C_residual << endl;
          Rcout<<"You can supress this waring by inrease the n.cell_type.residual to  " << C_all_current_raw - C_known << endl;
        }
        C_all_current_raw  = C_residual+C_known; // do not count the unknown cell type
      }
      
      
      
      MATTYPE U_all = find_gather_U(A,R_sample_all_cell_type_sum,C_all_current_raw,cell_type_min_cells,true,1e-8,false);
      C_all_current = U_all.n_rows;
      
      MATTYPE U_t = zeros<MATTYPE>(C_all_current, K);
      U_t.cols(all_clusters_select) = U_all;
      VECTYPE new_clusters_fraction = sum(U_t.cols(new_clusters), 1) / sum(U_all,1);
      float new_cell_type_new_clusters_fraction = 0.5;
      arma::uvec cell_types_index = sort_index(new_clusters_fraction, "descend");
      C_new_current = sum( new_clusters_fraction > new_cell_type_new_clusters_fraction );
      if(C_new_current<1){
        C_new_current = 1; 
      }
      
      MATTYPE U_new = U_all.rows(cell_types_index.head(C_new_current));
      
      
      // join cell types if centroid close too much, reduce the fake cell type
      // only use the main new cells to join cell types
      MATTYPE Phi_C_ = U_new * R_main_new_cell_type; 
      // should consider the rare cell type?
      MATTYPE Y_C_raw = Z_cos.cols(main_new_clusters_cells)*Phi_C_.t()+1e-9; // d x C_new_current
      MATTYPE Y_C,Cor_C;
      Y_C = arma::normalise(Y_C_raw,2,0);
      Cor_C = Y_C.t()* Y_C;
      Cor_C.diag().fill(-2.);
      while((Cor_C.max()>centroid_cor_threshold) & (C_new_current>=2)){
        uword im = Cor_C.index_max();
        uvec	sub = ind2sub( size(Cor_C), im);
        int c1 = sub.min(), c2 = sub.max();
        ROWVECTYPE u_combine = U_new.row(c1) + U_new.row(c2);
        VECTYPE y_raw_combine = Y_C_raw.col(c1) + Y_C_raw.col(c2);
        if(C_new_current >2){
          VECTYPE findex = ones<VECTYPE>(C_new_current);
          findex(c1) = 0;
          findex(c2) = 0;
          arma::uvec index = find(findex>0);
          U_new.head_rows(C_new_current - 2) = U_new.rows(index);
          Y_C_raw.head_cols(C_new_current - 2) = Y_C_raw.cols(index);
        }
        U_new.row(C_new_current - 2) = u_combine;
        Y_C_raw.col(C_new_current - 2) = y_raw_combine;
        C_new_current -=1;
        Y_C = arma::normalise(Y_C_raw.head_cols(C_new_current),2,0);
        Cor_C = Y_C.t()* Y_C;
        Cor_C.diag().fill(-2.);
      }
      U_new = U_new.head_rows(C_new_current);
      C_new = C_new_current;
      // layout C_known 1 C_new not need to ship the unknown column to the new positon
      C_current = C_known + C_new + 1;
      
      Phi_C.cols(main_new_clusters_cells).zeros();
      
      uvec cindex = regspace<arma::uvec>(C_known+1,C_current-1);
      Phi_C(cindex, main_new_clusters_cells) = 
        U_new * R_main_new_cell_type;
      Phi_C(cindex, main_new_clusters_cells) = normalise(Phi_C(cindex, main_new_clusters_cells), 1, 0);
      
      // // push to 1
      // arma::urowvec cell_type_position = index_max(Phi_C(cindex, main_new_clusters_cells),0);
      // ROWVECTYPE cell_type_prob = max(Phi_C(cindex, main_new_clusters_cells),0);
      // arma::uvec push_1_index = find(cell_type_prob > new_cell_type_prob_threshold);
      // cell_type_position = cell_type_position(push_1_index).t();
      // 
      // Phi_C(cindex, main_new_clusters_cells(push_1_index)).zeros();
      // for(int i =0; i<push_1_index.size();i++){
      //   Phi_C(cindex(cell_type_position(i)), main_new_clusters_cells(push_1_index(i)))=1.0;
      // }
      //Phi_C(cindex, main_new_clusters_cells) = normalise(Phi_C(cindex, main_new_clusters_cells), 1, 0);
      
      
      cells_type_known = unique(join_cols(cells_type_known_original, main_new_clusters_cells));
      VECTYPE unknown_ = ones<VECTYPE>(N);
      unknown_(cells_type_known).fill(0);
      cells_type_unknown = find(unknown_>0);
    }else{
      //Rcout <<"n_cells_U_new_select " <<n_cells_U_new_select<<"is too small. with new_cell_type_prob_threshold "<<new_cell_type_prob_threshold<<endl;
      C_new = 0;
      C_current = C_known+1;
      cells_type_known = cells_type_known_original;
      cells_type_unknown = cells_type_unknown_original;
    }
  }else{
    C_new = 0;
    C_current = C_known+1;
    cells_type_known = cells_type_known_original;
    cells_type_unknown = cells_type_unknown_original;
  }
  
  //Rcout<<" begin update_Phi_C infer "<<endl;
  // remove unknown column?
  Ru.elem(find(Ru< 1e-8)).fill(1e-8);
  Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) =
    (Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_known) % repmat(u_scale(cells_type_known).t(),C_current,1))
    * R.cols(cells_type_known).t() 
    % repmat(1./Ru.t(),C_current,1)
    * R.cols(cells_type_unknown);
    A_U = sqrt(alpha) * diagmat(1./sqrt(Ru)) * R.cols(cells_type_unknown);
    
    // VECTYPE eigval;
    // MATTYPE eigvec;
    // eig_sym(eigval, eigvec, A_U* A_U.t());
    // Rcout<<"eigvaule_max "<<max(eigval)<<endl;
    // float epilon_eigvalue = min(1.- max(eigval),1e-8);
    // Phi_C.cols(cells_type_unknown) = Phi_C.cols(cells_type_unknown) +
    //   Phi_C.cols(cells_type_unknown) * A_U.t() * 
    //   eigvec * arma::diagmat(1./(1. + epilon_eigvalue - eigval)) * eigvec.t() 
    //   * A_U;
    
    // VECTYPE Ru_U = R*u.t()*alpha;
    // float eigvaule_max_est = as_scalar(max(Ru_U/Ru));
    // float epilon_eigvalue = min(1.- eigvaule_max_est,1e-8);
    // Phi_C.cols(cells_type_unknown) = Phi_C.cols(cells_type_unknown) +
    //   Phi_C.cols(cells_type_unknown) * A_U.t() *  inv_sympd((1. + epilon_eigvalue) * arma::eye(K,K) - A_U* A_U.t()) * A_U;
    
    // using the solve bellow to avoid the singular problem eigvalue may near 0
    MATTYPE T_KC;
    T_KC = solve( eye(K,K) - A_U* A_U.t(), A_U*Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown).t());
    Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown)= Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) +
      T_KC.t() * A_U;
    Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown) = normalise(Phi_C(regspace<arma::uvec>(0,C_current-1),cells_type_unknown),1,0);
}

RCPP_MODULE(dincta_classification_module) {
  class_<dincta_classification>("dincta_classification")
  .constructor<int>()
  
  .field("Z_corr", &dincta_classification::Z_corr)
  .field("Z_orig", &dincta_classification::Z_orig)  
  .field("Z_cos", &dincta_classification::Z_cos)  
  .field("R", &dincta_classification::R)  
  .field("Y", &dincta_classification::Y)  
  .field("Phi_C", &dincta_classification::Phi_C) 
  .field("u", &dincta_classification::u)
  .field("objective_kmeans", &dincta_classification::objective_kmeans)
  .field("objective_kmeans_dist", &dincta_classification::objective_kmeans_dist)
  .field("objective_kmeans_entropy", &dincta_classification::objective_kmeans_entropy)
  .field("objective_kmeans_kl_cell_type_loss", &dincta_classification::objective_kmeans_kl_cell_type) 
  .field("centroid_max_cors", &dincta_classification::centroid_max_cors)
  .field("cell_type_min_cells_vector", &dincta_classification::cell_type_min_cells_vector)
  .field("dist_mat", &dincta_classification::dist_mat)
  .field("ran_setup", &dincta_classification::ran_setup)
  .field("ran_init", &dincta_classification::ran_init)
  
  
  .field("N", &dincta_classification::N)
  .field("K", &dincta_classification::K)
  .field("C", &dincta_classification::C)
  .field("C_known", &dincta_classification::C_known)
  .field("C_new", &dincta_classification::C_new)
  .field("d", &dincta_classification::d)
  .field("M", &dincta_classification::M)
  .field("max_iter_kmeans", &dincta_classification::max_iter_kmeans)
  //.field("frequency_update_Phi_C", &dincta_classification::frequency_update_Phi_C)
  
  .field("sigma", &dincta_classification::sigma)
  .field("sigma_entropy", &dincta_classification::sigma_entropy)
  .field("sigma_cell_type", &dincta_classification::sigma_cell_type)
  .field("mu", &dincta_classification::mu)
  .field("u_original", &dincta_classification::u_original)
  .field("O", &dincta_classification::O) 
  .field("update_order", &dincta_classification::update_order)    
  .field("cells_update", &dincta_classification::cells_update)    
  .field("kmeans_rounds", &dincta_classification::kmeans_rounds)    
  .field("update_Phi_C_rounds", &dincta_classification::update_Phi_C_rounds) 
  .field("epsilon_kmeans", &dincta_classification::epsilon_kmeans)    
  .field("epsilon_dincta", &dincta_classification::epsilon_dincta)
  
  //.method("init_cluster", &dincta_classification::init_cluster_cpp)
    .method("check_convergence", &dincta_classification::check_convergence)
    .method("setup", &dincta_classification::setup)
    .method("compute_objective", &dincta_classification::compute_objective)
    .method("update_R", &dincta_classification::update_R)
    .method("init_cluster_cpp", &dincta_classification::init_cluster_cpp)
    .method("init_refine_Phi_C_and_check_u_convergence_cpp", 
  &dincta_classification::init_refine_Phi_C_and_check_u_convergence_cpp)
  .method("cluster_cpp", &dincta_classification::cluster_cpp)
  
  ;
}




