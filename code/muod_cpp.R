source("cpp_functions/cluster_invoke_cpp_functions.R")

# This is the code to run Fast-MUOD and Semifast-MUOD as described in the paper
# Detecting and Classifying Outliers in Big Functional Data. 
# Main funtion to run is the get_outliers function. The parameters are:
# data: d by n matrix where n is the sample size and d is the evaluation points
# method: can take "rcpp" for MUOD, "fast" for Fast-MUOD with point-wise median,
#         "fast_l1" for Fast-MUOD with L1 median, and "semifast" for Semifast-MUOD
# cut_method: method to cut in the muod indices. Can take "boxplot" for boxplot which is 
#             used in paper, "adjboxplot" for skewness adjusted boxplot of Hubert et al. (2015)
#             "tangent" for the tangent method used in Azcorra et al. (2018), "carling_box" 
#             for the carling boxplot
# plot_cutoff: plot the sorted muod indices together with the selected cutoff
# sample_prop: sample proportion to use if "semifast" is selected for method. Should be a 
#              value between 0 and 1. Default is 0.5
# n_core: number of core to use for parallel processing. Advisable for method = "rcpp".
#          requires the parallel package.

get_outliers <- function(data, method = c("rcpp", "fast", "fast_l1", "semifast"),
                         cut_method = c("boxplot", "adjboxplot", "tangent", "log_transform",
                                        "carling_box", "carling_combo"),
                         plot_cutoff = F, sample_prop = 0.5, n_core = 3) {
  
  #setup params
  n <- ncol(data)
  options("mc.cores" = n_core)
  n_groups =  n/n_core # n/getOption("mc.cores", 1)
  benchmark = c(1, 0, 1)
  
  method <- match.arg(method)
  cut_method <- match.arg(cut_method)
  
  # import cpp implementation
  data <- sanitizeData(data) 
  
  # get muod indices
  indices <- get_muod_indices(data = data, n = n, n_groups = n_groups,
                              method = method, benchmark = benchmark,
                              sample_prop = sample_prop)
  
  # compute outliers cutoff
  get_outliers_from_indices(indices, cut_method, plot_cutoff, data,
                         c_value, n_sim, n_search, fdr, method)
}


get_muod_indices <- function(data, n = ncol(data), n_groups = n/getOption("mc.cores", 1),
                             method = c("rcpp", "fast", "semifast", "fast_l1"),
                             benchmark=c(1, 0, 1), sample_prop = 0.5, seed) {
  
  method <- match.arg(method)
  
  # initialize FORK cluster
  # FORK cluster necessary so that clusters have access to compiled cpp functions
  # without the need for explicit export. FORK doesn't work on windows.
  # this is a workaround till we build a package.
  par_cluster <- makeCluster(getOption("mc.cores", 1), type = "FORK") 
  
  # needed if using psock clusters
  #clusterEvalQ(par_cluster, source("../cpp_functions/cluster_invoke_cpp_functions.R"))
  #clusterExport(par_cluster, varlist = list("cor_cov_blockwise", "cov_armadillo",
  #                                          "cov_armadillo_mclsm_sf", "cov_armadillo_mclsm_normal"))
  
  if(method == "rcpp" ){
    pre_indices <- compute_muod_indices_rcpp(data, n, n_groups, par_cluster)
  } else if(method == "fast") {
    pre_indices <- compute_muod_indices_fast(data, n, n_groups, par_cluster)
  } else if(method == "semifast"){
    pre_indices <- compute_muod_indices_semifast(data, n, n_groups, par_cluster,
                                                 sample_prop, seed = seed)
  } else if(method == "fast_l1"){
    pre_indices <- compute_muod_indices_fast_l1(data, n, n_groups, par_cluster)
  } else {
    stop("not a valid MUOD method")
  }
  
  # apply benchmark
  abs(as.data.frame(pre_indices - matrix(benchmark, nrow(pre_indices), 3, byrow = T)))
}

###Rcpp Luisfo#####
compute_muod_indices_rcpp <- function(data, n = ncol(data),
                                      n_groups = n/getOption("mc.cores", 1),
                                      par_cluster){
  # compute the column splits/partition for parallel processing
  
  splits <- parallel:::splitList(1:n, max(1, as.integer(n/n_groups)))
  # compute auxiliar support data
  data_means <- colMeans(data, na.rm = T)
  data_vars <- apply(data, 2, var, na.rm = T)
  data_sds <- apply(data, 2, sd, na.rm = T)
  data2 <- t(t(data) - data_means) #pre computed mean-distance data
  # compute indices 
  #sourceCpp("cor_cov_blockwise.cpp") #in case
  vectors <- do.call(rbind, parLapply(par_cluster, splits, mean_cor_lsm_rcpp,
                                      data2, data_means, data_vars, data_sds))
  stopCluster(par_cluster)
  vectors <- data.frame(vectors)
  colnames(vectors) <- c("shape", "magnitude", "amplitude")
  vectors
}
mean_cor_lsm_rcpp <- function(i, mtx2, means, vars, sds){
  pre_outl <- cor_cov_blockwise(mtx2, means, vars, sds, i[1], i[length(i)])
  pre_outl[,1] <- pre_outl[, 1]/sds[i]
  pre_outl[,2] <- means[i] - pre_outl[,2]
  pre_outl
}

#######fast#############
compute_muod_indices_fast <- function(data, n = ncol(data),
                                      n_groups = n/getOption("mc.cores", 1),
                                      par_cluster){
  # compute reference observation (median of each var)
  data_ref <- apply_median_cpp(data, 1) 
  data_ref_mean <- mean(data_ref, na.rm = T)
  data_ref_var <- var(data_ref, na.rm = T)
  data_ref_sds <- sd(data_ref, na.rm = T)
  data_ref_dev <- data_ref - data_ref_mean
  # compute the column splits/partition for parallel processing
  splits <- parallel:::splitList(1:n, max(1, as.integer(n/n_groups)))
  # auxiliar vars
  data_means <- colMeans(data, na.rm = T)
  data_sds <- apply_sd_cpp(data, 2)
  data_dev <- t(t(data) - data_means)
  # compute the outlier values
  # compute Outliers
  vectors <- do.call(rbind, parLapply(par_cluster, splits, mean_cor_lsm_fast_rcpp, data_dev, data_ref_dev,
                                      data_ref_mean, data_ref_var, data_ref_sds,
                                      data_means, data_sds))
  stopCluster(par_cluster)
  vectors <- data.frame(vectors)
  colnames(vectors) <- c("shape", "magnitude", "amplitude")
  vectors
}

compute_muod_indices_fast_l1 <- function(data, n = ncol(data),
                                         n_groups = n/getOption("mc.cores", 1),
                                         par_cluster){
  
  # compute reference observation (using the L1median of data)
  data_ref <- pcaPP::l1median_NLM(t(data), maxit = 500)$par
  #data_ref <- apply_median_cpp(data, 1) 
  data_ref_mean <- mean(data_ref, na.rm = T)
  data_ref_var <- var(data_ref, na.rm = T)
  data_ref_sds <- sd(data_ref, na.rm = T)
  data_ref_dev <- data_ref - data_ref_mean
  # compute the column splits/partition for parallel processing
  splits <- parallel:::splitList(1:n, max(1, as.integer(n/n_groups)))
  # auxiliar vars
  data_means <- colMeans(data, na.rm = T)
  data_sds <- apply_sd_cpp(data, 2)
  data_dev <- t(t(data) - data_means)
  # compute the outlier values
  # compute Outliers
  vectors <- do.call(rbind, parLapply(par_cluster, splits, mean_cor_lsm_fast_rcpp, data_dev, data_ref_dev,
                                      data_ref_mean, data_ref_var, data_ref_sds, data_means, data_sds))
  stopCluster(par_cluster)
  vectors <- data.frame(vectors)
  colnames(vectors) <- c("shape", "magnitude", "amplitude")
  vectors
}

mean_cor_lsm_fast_rcpp <- function(i, mtx, mtx_ref, ref_mean,
                                   ref_var, ref_sds, means, sds){
  pre_outl <- corCovAggFast(mtx = mtx, mtxref = mtx_ref,
                          refmean = ref_mean, refvar = ref_var,
                          refsds = ref_sds, colstart = i[1],
                          colend = i[length(i)])
  
  pre_outl[,1] <- pre_outl[, 1]/sds[i]
  pre_outl[,2] <- means[i] - pre_outl[,2]
  pre_outl
}


  
  
#######semi fast #########

compute_muod_indices_semifast <- function(data, n = ncol(data),
                                          n_groups = n/getOption("mc.cores", 1),
                                          par_cluster, sample_prop = 0.5, seed){
  # check that prop is not more than 0.5
  
  # check if seed is miisng 
  if(missing(seed)){
    data_ref_index <- sample(1:n, size = ceiling(sample_prop * n)) 
    
  } else{
    set.seed(seed = seed)
    data_ref_index <- sample(1:n, size = ceiling(sample_prop * n)) 
    
  }
  # sample reference observations 
  
  data_ref <- data[, data_ref_index]
  
  
  # compute aux support data for ref data 
  data_ref_mean <- colMeans(data_ref, na.rm = T)
  data_ref_var <- apply_var_cpp(data_ref, 2)
  data_ref_sds <- apply_sd_cpp(data_ref, 2)
  data_ref_dev <- t(t(data_ref) - data_ref_mean)
  
  # compute the column splits/partition for parallel processing
  splits <- parallel:::splitList(1:n, max(1, as.integer(n/n_groups)))
  # auxiliar vars
  data_means <- colMeans(data, na.rm = T)
  data_sds <- apply_sd_cpp(data, 2)
  data_dev <- t(t(data) - data_means)
  # compute the outlier values
  
  # check that if ref data is jus one observation, need to readjust
  # dimension like in fast moud before sending to arma cov
  vectors <- do.call(rbind, parLapply(par_cluster, splits, mean_cor_lsm_semifast_rcpp_arma, data_dev, data_ref_dev,
                                      data_ref_mean, data_ref_var, data_ref_sds, data_means, data_sds))
  stopCluster(par_cluster)
  vectors <- data.frame(vectors)
  colnames(vectors) <- c("shape", "magnitude", "amplitude")
  vectors
}

mean_cor_lsm_semifast_rcpp_arma <- function(i, mtx, mtx_ref, ref_mean, ref_var,
                                            ref_sds, means, sds){
  pre_outl <- corCovAggSf(mtx = mtx, mtxref = mtx_ref,
                          refmean = ref_mean, refvar = ref_var,
                          refsds = ref_sds, colstart = i[1],
                          colend = i[length(i)])
  pre_outl[,1] <- pre_outl[, 1]/sds[i]
  pre_outl[,2] <- means[i] - pre_outl[,2]
  pre_outl
}













################################# cutting methods #################
# Cutting method
get_outliers_from_indices <- function(indices,
                                      cut_method = c("boxplot", "adjboxplot",
                                                     "bootstrap", "tangent",
                                                     "log_transform", "carling_box",
                                                     "carling_combo"),
                                      plot_cutoff = FALSE, data, c_value,
                                      n_sim, n_search, fdr, method = "") {
  cut_method <- match.arg(cut_method)
  # initialize clusters to use
  outl_type <- c("shape", "amplitude", "magnitude")
  sapply(outl_type, function(outl_name, plot_cutoff = F) { # think about plot param here
    # sort the curve
    metric <- sort(indices[, outl_name])
    # compute elbow point
    cutoff <- get_outlier_cutoff(metric, cut_method = cut_method,
                                 curve_name = outl_name,
                                 plot_cutoff = plot_cutoff,
                                 method = method)
    # filter outliers
    outl <- which(indices[, outl_name] > cutoff)
    
    outl[order(indices[outl, outl_name], decreasing = T)]
  }, plot_cutoff = plot_cutoff, simplify = F, USE.NAMES = T)
}


#######legacy cutoffs##################

get_outlier_cutoff <- function(curve, cut_method = c("boxplot", "tangent",
                                                 "adjboxplot", "log_transform",
                                                 "carling_box", "carling_combo"),
                               curve_name = "", plot_cutoff = FALSE, method = ""){
  cut_method <- match.arg(cut_method)
  if(cut_method == "adjboxplot"){
    box_stats <- robustbase::adjbox(x = curve, plot = F) 
    which_best <- which(curve == box_stats$stats[5, 1])
  } else if(cut_method == "boxplot"){
    box_stats <- boxplot.stats(curve)
    which_best <- which( curve == box_stats$stats[5])
  } else if(cut_method == "carling_box") {
    box_stats <- WRS::idealf(curve)
    nn <- length(curve)
    k = ((17.63*nn) - 23.64)/((7.74*nn) - 3.71)
    which_best <- median(curve) + (k*(box_stats$qu - box_stats$ql))
    return(which_best)
  }else if(cut_method == "carling_combo") {
    if(curve_name == "amplitude"){
      box_stats <- boxplot.stats(curve)
      which_best <- which( curve == box_stats$stats[5])
    }else{
      box_stats <- WRS::idealf(curve)
      nn <- length(curve)
      k = ((17.63*nn) - 23.64)/((7.74*nn) - 3.71)
      which_best <- median(curve) + (k*(box_stats$qu - box_stats$ql))
      return(which_best)
    }

  }else if(cut_method == "tangent"){
    which_best <- compute_tangent_cutoff(curve, plot = plot_cutoff)
  } else if(cut_method == "log_transform"){
    which_best <- compute_log_cutoff(curve, method = method)
  }else {
    stop("not a valid cutting method")
  }
  cutoff <- curve[which_best]
  if(plot_cutoff){
    
    plot(curve, ylab = paste0("sorted ", curve_name, " indices"),
         xlab = "n", main = paste0(curve_name, " index cut @ ", round(cutoff, 4)))
    abline(h = cutoff)
    
    print(ggplot(data = data.frame(metric = curve) ,
                 aes(y = sort(metric), x = 1:length(metric))) +
            geom_point(size =  .2) + 
            geom_hline(mapping = aes(yintercept = cutoff),
                       colour = "red") + 
            ggtitle(paste0("Indices of ", curve_name)) +
            ylab(paste0("sorted ", curve_name, " indices")) +
            xlab("n") +
            theme_bw()
    )
  }
  
  
  
  cutoff
}


## log_transform cutoff##

compute_log_cutoff <- function(metric, method = ""){
  quan <- qnorm(.995)
  if(method == "fast" | method == "fast_l1"){
    quan <- qnorm(.99)
    #cat("using 98.5 % because method is fastMUOD")
    }
  cpoint <- exp(median(log(metric)) + mad(log(metric), constant = 1)* quan)
  #fixplottin
  # if(plott){
  #  print(ggplot(data = data.frame(metric = metric) ,
  #          aes(y = sort(metric), x = 1:length(metric))) +
  #     geom_point(size =  .2) + 
  #     geom_hline(mapping = aes(yintercept = cpoint),
  #                colour = "red") + 
  #     geom_hline(mapping = aes(yintercept = metric[max(which(metric < cpoint))]),
  #                 colour = "blue") + 
  #     ggtitle(paste0("Indices of ", metric_name))
  #    )
  # }
  
  return(max(which(metric < cpoint)))
}

####tangent cutoff############
compute_tangent_cutoff <- function(metric, plot = F){
  cpoint <- find_tangent_X_intercept(seq(0, 1, length = length(metric)), metric, plot = plot)
  ceiling(cpoint * length(metric))
}
find_tangent_X_intercept <- function(x, y, which_x = which.max(diff(y)) + 1, plot = F){
  spl <- smooth.spline(y ~ x)
  newx_0 <- x[which.min(diff(diff(y)))]
  newx_1 <- mean(x[c(which_x - 1, which_x)], na.rm = T)
  pred0 <- predict(spl, x = newx_0, deriv = 0)
  pred1 <- predict(spl, x = newx_1, deriv = 1)
  
  #slope correction
  pred1$y <- max(pred1$y, 1.5)
  
  y_intercept <- pred0$y - (pred1$y*newx_0)
  x_intercept <- -y_intercept/pred1$y
  
  if(plot) {
    plot(x, y, type = "l", ylim = c(0, max(y)))
    abline(h = 0, col = 8)
    lines(spl, col = 2) # spline
    points(pred0, col = 2, pch = 19) # point to predict tangent 
    lines(x, (y_intercept + pred1$y * x), col = 3) # tangent (1st deriv. of spline at newx)
    points(x_intercept, 0, col = 3, pch = 19) # x intercept
  }
  x_intercept
}


