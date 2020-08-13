# functions to source and Rcpp codes and pacakges
library(Rcpp)
library(parallel)
# import cpp functions
sourceCpp("cpp_functions/apply_mean_cpp.cpp")
sourceCpp("cpp_functions/apply_median_cpp.cpp")
sourceCpp("cpp_functions/apply_sd_cpp.cpp")
sourceCpp("cpp_functions/apply_var_cpp.cpp")
sourceCpp('cpp_functions/sanitize_data.cpp')
##list rcpp files to source. 
sourceCpp("cpp_functions/cor_cov_blockwise.cpp") # mean_cor_lsm_rcpp by luisfo
sourceCpp("cpp_functions/cov_agg_mclsm_sf.cpp")  # mean_cor_lsm_semifast_rcpp_arma
sourceCpp("cpp_functions/cov_agg_mclsm_fast.cpp") # mean_cor_lsm_fast_rcpp


# dev
#sourceCpp("cpp_functions/cov_armadillo.cpp") # mean_cor_lsm_fast_rcpp2, testing
#sourceCpp("cpp_functions/cov_armadillo_mclsm_sf.cpp")  # mean_cor_lsm_semifast_rcpp_arma2, testing
#sourceCpp("cpp_functions/cov_armadillo_mclsm.cpp") # mean_cor_lsm_fast_rcpp_arma
#sourceCpp("cpp_functions/cov_armadillo_mclsm_normal.cpp") # mean_cor_lsm_normal_rcpp_arma