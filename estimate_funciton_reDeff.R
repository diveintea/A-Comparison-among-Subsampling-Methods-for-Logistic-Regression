################ compute reDeff ###################

## Input
# 1. x    -- predictors with labels (current subsample data)
# 2. x_c  -- predictors without labels in candicate set
# 3. beta -- parameters of Logistic regression
## Output
# reDeff for each data points in x_c

compute_reDeff = function(x, x_c, beta){
  prob0 = 1 - 1/(1+exp(x %*% t(t(beta)) ) )
  W0 = diag(as.vector(prob0 * (1-prob0)))
  M0 = t(x) %*% W0 %*% x
  d = rep(NA, nrow(x_c))
  for(i in 1:nrow(x_c)){
    x_new = t(x_c[i,])
    prob = as.vector(1 - 1/(1+exp(x_new %*% t(t(beta)) ) ) )
    w = prob * (1-prob)
    M =  w * t(x_new) %*% x_new
    d[i] = det(M0+M)
  }
  return(d)
}

################ do predict ###################
predict_mine = function(x_test, y_test, parameter, threshold = 0.5){
  z_pre = x_test %*% t(t(parameter))
  p_pre = 1 - 1/(1+exp(z_pre))
  y_pre = rep(0,length(p_pre))
  y_pre[p_pre > threshold] = 1
  return(list("prediction" = y_pre, "acc" = mean(y_test == y_pre)))
}



#################### reDeff subsampling####################

## Input
# 1. X      -- predictors
# 2. y      -- label, response
# 3. sub_r  -- subsample index in step 1 with size r0
#             (would be extended in function until length(sub_r)=r0+r1)
# 4. r1_set -- subsample size in step 2 (could be a vector)
# 5. q      -- target location
# 6. h      -- num. of data points near q
# 7. beta_tilde_0 -- parameters estimated on x[sub_r,-1]
## Output
# 1. beta_tilde -- the estimates on r0+r1 subsamples
# 2. sub_r1 -- indices of subsample in step2
# 3. sub_r0 -- indices of subsample in step1
# 4. sub_r1_situation -- record the lower and upper bound of candidate sets and selected points in each loop 
# 5. model_list -- a list of models. you can use ablities of glm(), but it will cost huge resource to store them.  

reDeff_subsample = function(x, y, sub_r0, r1_set, q, h, beta_tilde_0){
  ## Unifrom subsample and estimate beta_tilde_0
  n = nrow(x)
  beta_tilde_t = beta_tilde_0
  ## start reDeff subsample 
  sub_r = sub_r0
  sub_r1 = c()
  sub_r1_situation = c()
  ## initialize return variable
  beta_final = c()
  # model_list = list()
  
  for(t in 1:max(r1_set)) {#max(r1_set)
    ## find the set which includes the first h's data points near q
    ## get the indices of unlabel datas
    sub_r_c = c(1:n)[-sub_r]
    ## calculate response probabilities of unlabel data based on current beta_tilde_t
    prob = 1 - 1/(1+exp(x %*% beta_tilde_t))
    ## get the set 
    candidate = c()
    candidate_range = c()
    for(quantile in q){
      ## compute distance
      d = cbind(abs(prob[sub_r_c] - quantile), sub_r_c)
      # sort by distance
      d = d[order(d[,1]), ]
      ## take h/length(q)'s points into candidate for each q
      prob_this_candidate = prob[d[1:round(h/length(q)),2]]
      candidate_range = c(candidate_range, 
                          min(prob_this_candidate),
                          max(prob_this_candidate))
      candidate = c(candidate, d[1:round(h/length(q)), 2])
    }
    ## compute relative D-efficiency on candidate set
    reDeff = compute_reDeff(x[sub_r,], x[candidate,], beta_tilde_t)
    new_point = candidate[which.max(reDeff)]
    sub_r1 = c(sub_r1, new_point)
    sub_r = c(sub_r0, sub_r1)
    sub_r1_situation = cbind(sub_r1_situation, c(prob[new_point],candidate_range))
    ## re-estimate beta on the new sub_r
    # currenttime = proc.time()
    lr0 = glm(y[sub_r]~x[sub_r,-1],family = binomial(link = logit))
    beta_tilde_t = lr0$coefficients
    # cat(r1, " re-estimate beta_tilde_t, spend ", proc.time()-currenttime, "\n")
    if(t %in% r1_set){
      beta_final = rbind(beta_final, beta_tilde_t)
      model_list[[which(t == r1_set)]] = lr0
    }
  }
  return(list("beta" = beta_final, "sub_r1" = sub_r1, "sub_r1_situation" = sub_r1_situation, 
              "sub_r0" = sub_r0))
              # ,"model_list" = model_list))
}


############# verify the results in D-optimal Design of logistic regression #########

#################### compute D-optimal quantile ###################

## Input
# 1. beta0 -- the intercept of logistic regressoin
# 2. k     -- num. of predictors (except intercept term)
## Output
# two points which are D-optimal for logistic model
# (use logistic function F_L(h) to convert to prob.)
# 用很小的beta0解(b)即可

find_h = function(beta0, k){
  if(beta0 >= -1.5434){ # 不用
    fun = function(h){
      lamda_dot = (1-exp(h))/(1+exp(h))
      return(2 + lamda_dot*(h-beta0))
    }
    #curve(fun(x), beta0 , beta0+5)
    #abline(a = 0, b = 0, col = 'red')
    h_a = uniroot(fun, c(beta0 , beta0+5))$root
    
    return( c(beta0, h_a) )
  }else{
    fun = function(h){
      lamda_dot = (1-exp(h))/(1+exp(h))
      gamma_h = sqrt((2*k*beta0*h)^2+(h^2-beta0^2)^2)
      part1 = (h^2)*(2*k-1+k*lamda_dot*h)
      part2 = (beta0^2)*(1+k*lamda_dot*h)
      part3 = gamma_h*(1+lamda_dot*h)
      
      return(part1+part2+part3)
    }
    #curve(fun(x), 0 , -beta0)
    #abline(a = 0, b = 0, col = 'red')
    h_b = uniroot(fun, c(0,-beta0))$root
    
    return( c(-h_b, h_b) )
  }
}
################plot beta0 vs. p by dimention#######################
k = 1
h = c()
t = c()
for(b0 in seq(from=-5,to=1,by=0.001)){
  result = find_h(b0, k)
  h = rbind(h,result)
  t = rbind(t, result - b0)
}
plot(c(seq(from=-5,to=1,by=0.001)),t[,1], type = 'l', ylim = c(0,6))
lines(c(seq(from=-5,to=1,by=0.001)),t[,2],type = 'l')

1-1/(1+exp(h))[1,]

k_set = c(1,2,3,4,5,6,7,8)
b0_points = seq(from=-7,to=-1.534,by=0.001)
h = matrix(NA, ncol = length(b0_points), nrow = length(k_set))
for(i in 1:length(k_set)){
  for(j in 1:length(b0_points)){
    result = find_h(b0_points[j], k_set[i])
    h[i,j] = result[2]
  }
}

plot(b0_points,plogis(h[1,]),type = 'l', ylim = c(0,1))
lines(b0_points,plogis(-h[1,]))
for(i in 2:9){
  lines(b0_points,plogis(h[i,]), col=i)
  lines(b0_points,plogis(-h[i,]), col=i)
}
