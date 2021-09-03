library(wordspace)

################# compute SSP without labels ##################
## Input
# 1. prob_log_pilot -- response prob. of all data based on beta_pilot
# 2. X              -- predictors of all data
# 3. case           -- choice which case do you want to compute, mEVc or mEMSE
## Output
# SSP without labels
compute_SSP = function(prob_logi_pilot,X, case = "mEVc"){
  weight_pilot = prob_logi_pilot*(1-prob_logi_pilot)
  if(case == "mEVc"){
    # mEVc
    norm_xi = rowNorms(X)
    pi_mVc = (sqrt(weight_pilot) * norm_xi)/sum( sqrt(weight_pilot) * norm_xi)
    return(pi_mVc)
  }
  else{
    # mEMSE
    n = nrow(X)
    Mx_pilot = (1/n) * t(X) %*% (as.vector(weight_pilot) * X)
    Mx_pilot_inv = solve(Mx_pilot)
    norm_Mx_inv_xi = colNorms(Mx_pilot_inv %*% t(X))
    pi_mMSE = ( sqrt(weight_pilot) * norm_Mx_inv_xi)/sum( sqrt(weight_pilot) * norm_Mx_inv_xi)
    return(pi_mMSE)
  }
}


###################### OSMAC by bach  ########################
## Input
# 1. x -- predictors matrix (sample size) X (1+d)
# 2. y -- response (labels)
# 3. sub_r0 -- indices of subsample in step1 of size r0
# 4. r1_set -- subsample size in step2 (could be a vector)
# 5. beta_tilde_0 -- estimate on sub_r0
# 6. case -- indicate which case of SSP would be adopted (mEVc or mEMSE)
## Output
# 1. beta -- the final estimate on total subsample of size r0+r1 
#            (if r1 is a vector, beta would be a matrix (length(r1)) X (1+d))
# 2. sub_r1 -- indices of subsample in step2
# 3. sub_r0 -- indices of subsample in step1
# 4. model_list -- a list of models. you can use ablities of glm(), but it will cost huge resource to store them.  

OSMAC_subsample = function(x, y, sub_r0, r1_set, beta_tilde_0, case = 'mEVc'){
  ## Initailize or setup some variables
  n = nrow(x)
  beta_tilde_t = beta_tilde_0
  sub_r = sub_r0
  sub_r1 = c()
  beta_final = c()
  # model_list = list()
  
  ## start OSMAC subsample 
  for(t in 1:max(r1_set)) {
    ## compute SSP for the points which are not in sub_r
    ## get the indices of unlabel datas
    sub_r_c = c(1:n)[-sub_r]
    ## calculate response probabilities of unlabel data based on current beta_tilde_t
    prob = 1 - 1/(1+exp(x %*% beta_tilde_t))
    ## pi_opt (case = mEVc or mEMSE)
    pi_opt = compute_SSP(prob, x, case)
    pi_and_index = cbind(pi_opt,c(1:n))[sub_r_c,]
    new_point_index = pi_and_index[which.max(pi_and_index[,1]),2]
    
    sub_r1 = c(sub_r1, new_point_index)
    sub_r = c(sub_r0, sub_r1)

    ## re-estimate beta on the new sub_r
    lr0 = glm(y[sub_r]~x[sub_r,-1],family = binomial(link = logit))
    beta_tilde_t = lr0$coefficients
    ## check whether beta_tilde_t should be saved
    if(t %in% r1_set){
      beta_final = rbind(beta_final, beta_tilde_t)
      # model_list[[which(t == r1_set)]] = lr0
    }
  }
  return(list("beta" = beta_final, "sub_r1" = sub_r1 , "sub_r0" = sub_r0))
              # ,"model_list" = model_list))
}


#################### compute the MSE in Array ############################
## Input
# 1. beta_true    -- the object to compare with 
# 2. beta_3dArray -- the estimations in (experiments, predictors, subsample size)
## Output
# MSE from beta_3dArray to beta_true

get_MSE = function(beta_true, beta_3dArray){
  if(is.na(dim(beta_3dArray)[3])){
    mse_set = mean(apply((beta_3dArray - c(beta_true))^2, 1, sum))
  }else{
    num = dim(beta_3dArray)[3]  
    mse_set = rep(NA,num)
    for(i in 1:num){
      mse_set[i] = mean(apply((beta_3dArray[,,i] - c(beta_true))^2, 1, sum))
    }
  }
  return(mse_set)
}

#################### adjust SSP in finite sample ########################
## Input
# 1. r0       -- subsample size in step-1
# 2. r1       -- subsample size in step-2
# 3. pi_opt   -- optimal SSP a.k.a. final proportion of each samples
# 4. pi_step1 -- SSP in step-1 (usually is 1/n, the uniform SSP)
## Output
# SSP in step-2 such that the final proportion of each samples would approximate pi_opt
adjust_SSP = function(r0,r1,pi_opt,pi_step1){
  r = r0+r1
  frac = (1/r1)*(r*pi_opt - r0*pi_step1)
  frac[frac < 0] = 0
  pi_step2 = frac/ sum(frac)
  return(pi_step2)
}

################# estimate function ######################
## Input
# 1. y  -- response of subsamples
# 2. X  -- predictors of subsamples
# 3. pi -- the SSPs of each subsamples
#          a.k.a. reciprocal of weights in weighted MLE
# 4. beta_initial -- the initial point for optimizer
## Output
# the estimation of parameters in logistic regression by weighted MLE

weighted_logistic_nloptr = function(y,X,pi,beta_initial=NULL){
  # x = (intercept,x1,...,xd-1) 
  r = length(y)
  #設定目標函數為negative weighted log-likelihood
  nega_wei_log_likelihood = function(beta){
    beta = t(t(beta))
    prob_logi_beta = 1-1/(1+exp(X %*% beta))
    #prob_logi_beta[is.na(prob_logi_beta)] = 0.999999
    prob_logi_beta[prob_logi_beta > 0.999999] = 0.999999
    prob_logi_beta[prob_logi_beta < 0.000001] = 0.000001
    return(  -sum((1/(pi))*(y*log(prob_logi_beta)+(1-y)*log(1-prob_logi_beta)))     )
  }
  #設定目標函數的梯度函數
  nega_score = function(beta){
    beta = t(t(beta))
    prob_logi_beta = 1-1/(1+exp(X %*% beta))
    #prob_logi_beta[is.na(prob_logi_beta)] = 0.999999
    prob_logi_beta[prob_logi_beta > 0.999999] = 0.999999
    prob_logi_beta[prob_logi_beta < 0.000001] = 0.000001
    return(  -apply(as.vector((1/(pi))*(y-prob_logi_beta))*X,2,mean) )
  }
  
  if(is.null(beta_initial)){
    fullDataMatrix = as.data.frame(cbind(X,y))
    lr0 = glm(y~. ,data = fullDataMatrix[,-1] ,family = binomial(link = logit))
    beta_initial = lr0$coefficients
  }
  
  #設定演算法
  opt = list("algorithm"="NLOPT_LD_CCSAQ","xtol_rel"=1.0e-6) #NLOPT_LD_MMA
  #輸入起始點,目標函數,梯度,與演算法，開始估計beta
  res = nloptr( x0 = beta_initial,# rep(0.5,ncol(X)),
                eval_f = nega_wei_log_likelihood,
                eval_grad_f = nega_score,
                opts=opt)
  #iter_set = rbind(iter_set,res$iterations)
  #將此次求出的beta存入beta_set
  return(res$solution)
}
