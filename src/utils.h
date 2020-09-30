void cosine_normalize(MATTYPE& X, int margin, bool do_safe) {
  // to avoid Inf values, first divide by max 
  if (margin == 1) {
    for (unsigned r = 0; r < X.n_rows; r++) {
      if (do_safe)
        X.row(r) = X.row(r) / X.row(r).max(); 
      X.row(r) = X.row(r) / norm(X.row(r), 2);
    }     
  } else {
    for (unsigned c = 0; c < X.n_cols; c++) {
      if (do_safe)
        X.col(c) = X.col(c) / X.col(c).max(); 
      X.col(c) = X.col(c) / norm(X.col(c), 2);
    }    
  } 
}

//  negative entropy 
MATTYPE safe_entropy(const MATTYPE& X) {
  MATTYPE A = X % log(X);
  A.elem(find_nonfinite(A)).zeros();
  return(A);
}


// Overload pow to work on a MATTYPErix and vector
MATTYPE pow(MATTYPE A, const VECTYPE & T) {
  for (unsigned c = 0; c < A.n_cols; c++) {
    A.col(c) = pow(A.col(c), as_scalar(T.row(c)));
  }
  return(A);
}

MATTYPE find_gather_U(const MATTYPE &A, const VECTYPE K_sum, int C, float n_cell_type_min_cells, bool strict, float epsilon, bool verbose){
  int K = A.n_rows;
  MATTYPE U = eye(K,K);
  MATTYPE E = A - U;
  MATTYPE D = E.t()*E;
  VECTYPE err = D.diag();
  
  MATTYPE O = D;
  O.diag().fill(1e10);
  vector<float> loss;
  loss.push_back(as_scalar(accu(err)));

  // if(verbose){
  //   Rcout<< "number of cell types: "<<K<< ", err:"<< loss[loss.size() - 1] << endl;
  // }
  
  int C_find = K;
  float err_current = 1e10;
  
  
  VECTYPE err_join,findex;
  int c1 = -1, c2=-1,C_=K-1,n_max=0;
  uword select_c;
  uvec index,index_,join_index_,index_v;
  umat join_index;
  MATTYPE D_;
  ROWVECTYPE u_combine;
  for(C_ = K-1; C_ >= C; C_ = C_ -1){
    join_index_ = find(O == min(min(O)));
    join_index  = ind2sub( size(O), join_index_);
    n_max =join_index.n_cols;
    if(n_max>1){
      err_join = zeros<VECTYPE>(n_max) - 1e10;
      for(int ci = 0; ci < n_max; ci++ ){
        int irow = join_index(0,ci);
        int icol = join_index(1,ci);
        if(irow<icol){
          err_join[ci] = err[irow]+ err[icol];
        }
      }
      select_c =index_max(err_join);
      // c1 < c2
      c1 = join_index(0,select_c) >join_index(1,select_c)? join_index(1,select_c):join_index(0,select_c);
      c2 = join_index(0,select_c) >join_index(1,select_c)? join_index(0,select_c):join_index(1,select_c);
    }
    
    //combine D rows c1, c2, cols c1,c2  to D_ \in R^{C_ \times C_}
    findex = ones<VECTYPE>(C_+1);
    findex[c1] = 0;
    findex[c2] = 0;
    index_v =find(findex==1);
    
    D_ = zeros<MATTYPE>(C_,C_) + 1e10;
    if(C_>1){
      //check span(0,0) work
      D_(span(0, C_-2), span(0, C_-2)) = D(index_v,index_v);
      D_(span(0, C_-2), C_-1) = D(index_v,c1*ones<uvec>(1)) +  D(index_v,c2*ones<uvec>(1));
      D_(C_-1, span(0, C_-2)) = D(c1*ones<uvec>(1),index_v) +  D(c2*ones<uvec>(1),index_v);
    }
    D_(C_-1,C_-1) = D(c1,c1) + D(c2,c2) + D(c1,c2) + D(c2,c1);
    err = D_.diag();
    err_current = as_scalar(accu(err));
    if(((err_current >loss[loss.size() -1]) | (abs(err_current - loss[loss.size() -1])<epsilon))&(!strict)){
      break;
    }
    D = D_;
    O = D;
    O.diag().fill( 1e10 );
    loss.push_back(err_current);
   
    
    u_combine =  U.row(c1) + U.row(c2);
    if(C_>1){
      U.rows(0, C_-2) = U.rows(index_v);
    }
    U.row(C_-1) =  u_combine;
    C_find = C_;
    
    // if(verbose){
    //   // if(c1==c2){
    //   //   Rcout<<"err (c1,c2)"<<c1<<" "<<c2<<endl;
    //   // }
    //   // if(as_scalar(accu(U.rows(0, C_find-1) == 1)) != K){
    //   //   Rcout<<"err sum(U)"<<as_scalar(accu(U.rows(0, C_find-1) == 1)) <<endl;
    //   //   Rcout<<"(c1,c2)"<<c1<<" "<<c2<<endl;
    //   //   Rcout<< "D "<< D <<endl;
    //   //   Rcout<< "U "<< U.rows(0,C_find-1) <<endl;
    //   // }
    //     
    //   Rcout<< "number of cell types: "<<C_find<< ", err:"<< loss[loss.size() - 1]<< endl;
    // }
    // print(paste("U:"))
    // for(ic in 1:C_ ){
    //      print(paste(U[ic,], collapse = ', '))
    //  }
  }
  
  // join the cell type whose numuber of cells less than n_cell_type_min_cells
 
  // if(verbose){
  //   Rcout<< "number of cell types: "<<K<< ", err:"<< loss[loss.size() - 1] << endl;
  // }
  VECTYPE d,o;
  bool shrink = true;
  VECTYPE C_sum ;
  if(C_find<=1){
    shrink =false;
  }else{
    C_sum = U.rows(0,C_find-1) * K_sum;
    shrink = sum(C_sum < n_cell_type_min_cells)>0;
  }
  while(shrink){
    
    index = find(C_sum < n_cell_type_min_cells);
    d = D.diag();
    
    if(index.size()>1){
      // find the largeest error in D
      index_ = find(d(index) == max(d(index)));
      c1 = index(index_(0));
    }else{
      c1 = index(0);
    }
    //find which cell type c2 to join
    o = D.col(c1);
    o(c1) = 1e10;
    index = find(o == min(o));
    if(index.size()>1){
      // find the largeest error in D
      index_ = find(d(index) == max(d(index)));
      c2 = index(index_(0));
    }else{
      c2 = index(0);
    }
    C_ = C_find-1;
    findex = ones<VECTYPE>(C_+1);
    findex[c1] = 0;
    findex[c2] = 0;
    index_v =find(findex==1);
  
    D_ = zeros<MATTYPE>(C_,C_) + 1e10;
    if(C_>1){
      D_(span(0, C_-2), span(0, C_-2)) = D(index_v,index_v);
      D_(span(0, C_-2), C_-1) = D(index_v,c1*ones<uvec>(1)) +  D(index_v,c2*ones<uvec>(1));
      D_(C_-1, span(0, C_-2)) = D(c1*ones<uvec>(1),index_v) +  D(c2*ones<uvec>(1),index_v);
    }
    D_(C_-1,C_-1) = D(c1,c1) + D(c2,c2) + D(c1,c2) + D(c2,c1);
    err = D_.diag();
    err_current = as_scalar(accu(err));
    D = D_;
    loss.push_back(err_current);
    
    u_combine =  U.row(c1) + U.row(c2);
    if(C_>1){
      U.rows(0, C_-2) = U.rows(index_v);
    }
    U.row(C_-1) =  u_combine;
    C_find = C_;
    if(C_find<=1){
      shrink =false;
    }else{
      C_sum = U.rows(0,C_find-1) * K_sum;
      shrink = sum(C_sum < n_cell_type_min_cells)>0;
    }
    
    // if(verbose){
    //   // if(c1==c2){
    //   //   Rcout<<"err (c1,c2)"<<c1<<" "<<c2<<endl;
    //   // }
    //   // if(as_scalar(accu(U.rows(0, C_find-1) == 1)) != K){
    //   //   Rcout<<"err sum(U)"<<as_scalar(accu(U.rows(0, C_find-1) == 1)) <<endl;
    //   //   Rcout<<"(c1,c2)"<<c1<<" "<<c2<<endl;
    //   //   Rcout<< "D "<< D <<endl;
    //   //   Rcout<< "U "<< U.rows(0,C_find-1) <<endl;
    //   // }
    //   //Rcout<< "number of cell types: "<< C_find << ", err:"<< loss[loss.size() - 1] << endl;
    // }
  }
  if(verbose){
    Rcout << "Final, find number of cell types: "<<  C_find <<", err:" << loss[loss.size() -1] << endl; 
    if(C_find>C){
      Rcout <<"If you want to reduce to the number of cell types: "<< C << ", set the parameter strict=true." << endl;
    }
  }
  if(C_find>1){
    return U.rows(0,C_find-1);
  }else{
    return U.row(0);
  }
  
}


MATTYPE merge_R(const MATTYPE & R, float thresh = 0.8) {
  MATTYPE cor_res = cor(R.t());
  int K = R.n_rows;
  
  // define equivalence classes
  uvec equiv_classes = linspace<uvec>(0, K - 1, K);
  int new_idx;
  for (int i = 0; i < K - 1; i++) {
    for (int j = i + 1; j < K; j++) {
      if (cor_res(i, j) > thresh) {
        new_idx = min(equiv_classes(i), equiv_classes(j));
        equiv_classes(i) = new_idx;
        equiv_classes(j) = new_idx;
      }
    }
  }
  
  // sum over equivalence classes
  uvec uclasses = unique(equiv_classes);
  MATTYPE R_new = zeros<MATTYPE>(uclasses.n_elem, R.n_cols); 
  for (unsigned i = 0; i < R_new.n_rows; i++) {
    uvec idx = find(equiv_classes == uclasses(i)); 
    R_new.row(i) = sum(R.rows(idx), 0);
  }
  return R_new;  
}

// [[Rcpp::export]]
MATTYPE compute_Y(const MATTYPE& Z_cos, const MATTYPE& R) {
  return arma::normalise(Z_cos * R.t(), 2, 0);
}


// [[Rcpp::export]]
MATTYPE scaleRows_dgc(const VECTYPE& x, const VECTYPE& p, const VECTYPE& i, 
                      int ncol, int nrow, float thresh) {
  // (0) fill in non-zero elements
  MATTYPE res = arma::zeros<MATTYPE>(nrow, ncol);
  for (int c = 0; c < ncol; c++) {
    for (int j = p[c]; j < p[c + 1]; j++) {
      res(i[j], c) = x(j);
    }
  }
  
  // (1) compute means
  VECTYPE mean_vec = arma::zeros<VECTYPE>(nrow);
  for (int c = 0; c < ncol; c++) {
    for (int j = p[c]; j < p[c + 1]; j++) {
      mean_vec(i[j]) += x[j];
    }
  }
  mean_vec /= ncol;
  
  // (2) compute SDs
  VECTYPE sd_vec = arma::zeros<VECTYPE>(nrow);
  arma::uvec nz = arma::zeros<arma::uvec>(nrow);
  nz.fill(ncol);
  for (int c = 0; c < ncol; c++) {
    for (int j = p[c]; j < p[c + 1]; j++) {
      sd_vec(i[j]) += (x[j] - mean_vec(i[j])) * (x[j] - mean_vec(i[j])); // (x - mu)^2
      nz(i[j])--;
    }
  }
  
  // count for the zeros
  for (int r = 0; r < nrow; r++) {
    sd_vec(r) += nz(r) * mean_vec(r) * mean_vec(r);
  }
  
  sd_vec = arma::sqrt(sd_vec / (ncol - 1));
  
  // (3) scale values
  res.each_col() -= mean_vec;
  res.each_col() /= sd_vec;
  res.elem(find(res > thresh)).fill(thresh);
  res.elem(find(res < -thresh)).fill(-thresh);
  return res;
}