#include <TMB.hpp>

template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}



template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_VECTOR(catches);
  
  DATA_VECTOR(index);
  
  DATA_INTEGER(years);
  
  DATA_SCALAR(sigma_proc);
  
  DATA_SCALAR(sigma_catch);
  
  PARAMETER(log_r);
  
  PARAMETER(log_k);
  
  PARAMETER(log_index_q);
  
  PARAMETER(log_init_dep);
  
  PARAMETER(log_sigma_obs);
  
  PARAMETER_VECTOR(effish);
  
  PARAMETER_VECTOR(log_pros_error);
  
  int y1 = 0;

  y1 = years + 1;

  vector<Type> b(y1);

  vector<Type> b_bmsy(y1);
  
  vector<Type> u_umsy(y1);
  
  vector<Type> pros_error(years);
  
  vector<Type> catches_hat(years);

  vector<Type> index_hat(years);

  vector<Type> u(years);

  Type sigma_obs = exp(log_sigma_obs);

  Type r = exp(log_r);

  Type k = exp(log_k);

  Type index_q = exp(log_index_q);

  Type init_dep = exp(log_init_dep);
  // 
  Type nll = 0.0;

  Type pen = 0.0;

  Type eps = 0.001;

  b(0) = k * init_dep;
  // 
  for(int t=0; t<years; t++) {
    
    u(t) = 1.0 / (1.0+exp(-effish(t)));
    
    pros_error(t) = exp(log_pros_error(t) - pow(sigma_proc,2) / 2);
  
    index_hat(t) = b(t) * index_q;
    // 
    catches_hat(t) =  b(t) * u(t);
    // 
    b(t+1) =  (b(t) * (1 + r * (1 - b(t) / k)) - catches_hat(t)) * pros_error(t); // check this
    // 
    b(t+1) = posfun(b(t+1), eps, pen);
    
  }
  
  nll -= dnorm(log_pros_error, Type(0.0), sigma_proc, true).sum();
  
  nll -= dnorm(log(catches), log(catches_hat), sigma_catch, true).sum();
  
  nll -= dnorm(log(index), log(index_hat), sigma_obs, true).sum();
  
  vector<Type> log_index_hat = log(index_hat);
  
  b_bmsy = b / (k / 2);
    
  u_umsy = u / (r / 2);
    
  ADREPORT(index_hat);
  
  ADREPORT(log_index_hat);
  
  ADREPORT(b);
  
  ADREPORT(u)
    
  ADREPORT(b_bmsy);
  
  ADREPORT(u_umsy)
  
  ADREPORT(catches_hat);
  // 
  REPORT(u);
  
  return nll;
  
  /* Quick Reference
   ===============
   
   ** Macros to read data and declare parameters:
   
   _Template_Syntax_              _C++_type_                     _R_type_
   DATA_VECTOR(name)              vector<Type>                   vector
   DATA_MATRIX(name)              matrix<Type>                   matrix
   DATA_SCALAR(name)              Type                           numeric(1)
   DATA_INTEGER(name)             int                            integer(1)
   DATA_FACTOR(name)              vector<int>                    factor
   DATA_SPARSE_MATRIX(name)       Eigen::SparseMatrix<Type>      dgTMatrix
   DATA_ARRAY(name)               array<Type>                    array
   PARAMETER_MATRIX(name)         matrix<Type>                   matrix
   PARAMETER_VECTOR(name)         vector<Type>                   vector
   PARAMETER_ARRAY(name)          array<Type>                    array
   PARAMETER(name)                Type                           numeric(1)
   
   ** Macro to report intermediate expressions back to R:
   
   REPORT(x)
   ADREPORT(x)
   
   ** Basic constructors:
   
   vector<Type> v(n1);
   matrix<Type> m(n1,n2);
   array<Type> a(n1,n2,n3)
   
   ** Basic operations:
   
   v+v,v-v,v*v,v/v                Pointwise binary operations
   m*v                            Matrix-vector multiply
   a.col(i)                       R equivalent of a[,,i]
   a.col(i).col(j)                R equivalent of a[,j,i]
   a(i,j,k)                       R equivalent of a[i,j,k]
   exp(v)                         Pointwise math
   m(i,j)                         R equivalent of m[i,j]
   v.sum()                        R equivalent of sum(v)
   m.transpose()                  R equivalent of t(m)
   
   ** Distributions:
   
   Type dnbinom2(const Type &x, const Type &mu, const Type &var, int give_log=0)
   Type dpois(const Type &x, const Type &lambda, int give_log=0)
   Type dlgamma(Type y, Type shape, Type scale, int give_log=0)
   Type dnorm(Type x, Type mean, Type sd, int give_log=0)
   
   ** Parallel accumulator declaration (only methods "+=" and "-="):
   
   parallel_accumulator<Type> res(this);
   
   */
  
}

