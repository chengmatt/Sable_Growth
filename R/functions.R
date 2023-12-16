add_newton = function(n.newton, ad_model, mle_optim) {
  
  tryCatch(expr = for(i in 1:n.newton) {
    g = as.numeric(ad_model$gr(mle_optim$par))
    h = optimHess(mle_optim$par, fn = ad_model$fn, gr = ad_model$gr)
    mle_optim$par = mle_optim$par - solve(h,g)
    mle_optim$objective = ad_model$fn(mle_optim$par)
  }, error = function(e){e})
  
}

get_al_trans_matrix = function(age_bins, len_bins, mean_length, sd) {
  
  # Construct age length matrix
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))
  age_length = matrix(0.0, nrow = length(age_bins), ncol = length(len_bins))
  
  for(age_ndx in 1:length(age_bins)) {
    for(len_ndx in 1:length(len_bins)) {
      if (len_ndx == 1) {
        age_length[age_ndx, len_ndx] = pnorm(len_bins[2], mean_length[age_ndx], sd[age_ndx])
      } else if (len_ndx == length(len_bins)) {
        age_length[age_ndx, len_ndx] = 1 - pnorm(len_bins[length(len_bins)], mean_length[age_ndx], sd[age_ndx])
      } else {
        age_length[age_ndx, len_ndx] = pnorm(len_bins[len_ndx+1], mean_length[age_ndx], sd[age_ndx]) -  
                                           pnorm(len_bins[len_ndx], mean_length[age_ndx], sd[age_ndx])
      }
    }
  }
  # age_length = prop.table(age_length, margin = 1) # renomalize if they aren't

  return(age_length)
} # end function


prepare_data = function(data, obs_sd_mat, age_bins, growth_model_name, re_model_name, var_param, years, n_proj_years) {
  
  # growth models
  if(growth_model_name == 'length') growth_model = 1 # length-at-age
  if(growth_model_name == 'weight') growth_model = 2 # weight-at-age
  
  # random effects models
  if(re_model_name == 'iid') re_model = 0 # iid
  if(re_model_name == '2dar1') re_model = 2 # 2dar1
  if(re_model_name == '3dar1') re_model = 3 # 3dar1
  if(re_model_name == 'constant') re_model = 500 # 3dar1
  
  # get age year index for indexing in GMRF
  ay_Index <- as.matrix(expand.grid("age" = seq_len(length(unique(age_bins))), 
                                    "year" = seq_len(length(unique(years)) + n_proj_years ) ))
  
  # Set up parameters
  parameters = list(ln_X_inf = log(80), # asymptotic length or weight
                    ln_Lmin = log(40), # theoretical size at t0
                    ln_k = log(0.2), # growth rate
                    ln_beta = log(3.02), # allometric scaling
                    ln_obs_sigma2 = rep(log(4), 2), # variance for observations
                    ln_alpha = log(1e-5),
                    rho_a = 0, # correlation by age
                    rho_y = 0, # correlation by year
                    rho_c = 0, # correlation by cohort
                    ln_eps_at = array(data = 0, dim = c(length(unique(age_bins)), # process errors
                                                        length(unique(years)) + n_proj_years ) ),
                    ln_eps_sigma2 = log(1)) # variance of process errors
  
  # set up data
  data_list = list(obs_mat = as.matrix(data),
              obs_sd_mat = obs_sd_mat,
              ages = as.vector(sort(unique(age_bins))),
              ay_Index = ay_Index, var_param = var_param, 
              re_model = re_model, growth_model = growth_model)


# Mapping stuff -----------------------------------------------------------

  map = list(ln_beta = factor(NA)) # set up parameters to fit (fix beta at 3.02)
  
  if(growth_model_name == "length") {
    map$ln_alpha = factor(NA)
  } # length models
  
  if(growth_model_name == "weight") {
    map$ln_obs_sigma2 = factor(rep(NA, 2)) # map out observed sigma
    map$ln_X_inf = factor(NA) # Fixed at arbitrary value
  } # weight models
  
  if(re_model_name == "2dar1") {
    map$rho_c = factor(NA) # don't estimate cohort correlation for 2dar1
  }
  
  if(re_model_name == 'iid') {
    map$rho_a = factor(NA)
    map$rho_y = factor(NA)
    map$rho_c = factor(NA)
  } # if iid, don't estimate these
  
  if(re_model_name == "constant") { # constant (no time-variation)
    map$rho_a = factor(NA)
    map$rho_y = factor(NA)
    map$rho_c = factor(NA)
    map$ln_eps_sigma2 = factor(NA)
    map$ln_eps_at = factor(rep(NA, length(parameters$ln_eps_at)))
  }
  
  return(list(data = data_list, parameters = parameters, map = map))
} # end function