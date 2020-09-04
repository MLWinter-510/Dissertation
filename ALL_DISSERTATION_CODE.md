# Dissertation
All code to run the simulations in my Masters Dissertation

# Code to run Equal Randomisation Simulations

Equal_random = function(k, sample_size, block, no_treat, prob_vec){
  
  set.seed(k)
  
  treatment = vector()
  
  outcome = vector()
  
  no_s = matrix(0L, no_treat, 1)
  
  no_f = matrix(0L, no_treat, 1)
  
  results = matrix(0L, nrow = no_treat, ncol = 2)
  
  
  for (t in 0:((sample_size/block) - 1)) {
    
    treatment = sample(no_treat, sample_size, prob = rep(1/no_treat, no_treat), replace = TRUE) # equal probability of being assigned
    
    #  rep(1:no_treat, each = sample_size/no_treat) replace above with this if want exactly equal allocation
    
    for(p in 1:block){
    
    for (k in 1:no_treat){
    
    if(treatment[p + (t*block)] == k ){
    
    outcome[p + (t*block)] = sample(c(1,0), size = 1, replace = TRUE, prob = c(prob_vec[k], 1 - prob_vec[k]))
    
      } 
    
    }
    
    if(outcome[p + (t*block)] == 1){
    
    no_s[treatment[p + (t*block)]] = no_s[treatment[p + (t*block)]] + 1
    
    }else{
    
    no_f[treatment[p + (t*block)]] = no_f[treatment[p + (t*block)]] + 1} 
    
      } 
    
    }
    
    
  
  results[,1] = no_s[,1]
  
  results[,2] = no_f[,1]
  
  outcome = c(results[,1], results[,2])
  
  
  p.values = rep(NA, no_treat-1)
  
  for (i in 2:no_treat) {
  
  if(nrow(results) == 2){
  
  test.results = results
  
  }else{
  
  test.results = rbind(results[1,], results[i,])}
  
  test.res = fisher.test(test.results, alternative = "less")
  
  p.values[i-1] = test.res$p.value}
  
  output = list("success" = no_s, "failure" =  no_f, "results" = results, "p.values" = p.values)
  
  return(c(outcome, p.values))

}


# Code to run Thompson Sampling Simulations

Thompson_Sampling = function(k, sample_size, block, no_treat, prob_vec, prior_mat){

set.seed(k)

treatment = vector()

outcome = vector()

no_s = matrix(0L, no_treat, 1)

no_f = matrix(0L, no_treat, 1)

results = matrix(0L, nrow = no_treat, ncol = 2)


for(t in 0:((sample_size/block) - 1)){
    
    # Calculating the allocation probabilities for each treatment at the start of each block   
    
    if(no_treat == 2){
    
    theta_2 = rbeta(100, no_s[2] + prior_mat[2,1], no_f[2] + prior_mat[2,2])
    
    theta_1 = rbeta(100, no_s[1] + prior_mat[1,1], no_f[1] + prior_mat[1,2])
    
    theta = cbind(theta_2, theta_1)
    
    alloc_exp = sum(theta[,1] > theta[,2])/(nrow(theta))
    
    alloc_p = c(1-alloc_exp, alloc_exp)
    
    }else{
    
    theta_3 = rbeta(100, no_s[3] + prior_mat[3,1], no_f[3] + prior_mat[3,2])
    
    theta_2 = rbeta(100, no_s[2] + prior_mat[2,1], no_f[2] + prior_mat[2,2])
    
    theta_1 = rbeta(100, no_s[1] + prior_mat[1,1], no_f[1] + prior_mat[1,2])
    
    theta = cbind(theta_3, theta_2, theta_1)
    
    alloc_exp2 = sum((theta[,1] > theta[,2])&(theta[,1] > theta[,3]))/(nrow(theta)) # allocate treatment 2
    
    alloc_exp1 = sum((theta[,2] > theta[,1])&(theta[,2] > theta[,3]))/(nrow(theta)) # allocate teatment 1
    
    alloc_con = sum((theta[,3] > theta[,2])&(theta[,3] > theta[,1]))/(nrow(theta))  # allocate control  
    
    alloc_p = c(alloc_con, alloc_exp1, alloc_exp2)
    
    }
    
    
    for(p in 1:block){
    
    treatment[p + (t*block)] = sample(no_treat, size = 1, replace = TRUE, prob = alloc_p)
    
    for (k in 1:no_treat){
    
    if(treatment[p + (t*block)] == k ){
    
    outcome[p + (t*block)] = sample(c(1,0), size = 1, replace = TRUE, prob = c(prob_vec[k], 1-prob_vec[k]))}
    
    }
    
    if(outcome[p + (t*block)] == 1){
    
    no_s[treatment[p + (t*block)]] = no_s[treatment[p + (t*block)]] + 1
    
    }else{
    
    no_f[treatment[p + (t*block)]] = no_f[treatment[p + (t*block)]] + 1} 
    
      }
    }
    
    
  results[,1] = no_s[,1]
  
  results[,2] = no_f[,1]
  
  outcome = c(results[,1], results[,2])
  
  p.values = rep(NA, no_treat-1)
  
  for (i in 2:no_treat) {
  
  if(nrow(results) == 2){
  
  test.results = results
  
  }else{
  
  test.results = rbind(results[1,], results[i,])}
  
  test.res = fisher.test(test.results, alternative = "less")
  
  p.values[i-1] = test.res$p.value}
  
  output = list("success" = no_s, "failure" =  no_f, "results" = results, "p.values" = p.values)
  
  return(c(outcome, p.values))
}


# Function used to calculated the allocation probabilities at the start of each block based upon the Gittins Index

Gittins_mat = read.table("Gittins_mat.txt")

Input_alloc = function(info_mat, block, noRuns, no_cat, no_treat){

covariatedata = round(runif(noRuns, 1, no_cat))

selected = matrix(0L, noRuns, block + 1)

z = matrix(NA, noRuns, block + 1)


for (j in 1:noRuns){

s = info_mat[,1] + 1

f = info_mat[,2] + 1

for(t in 0:(block - 1)){

index = c()

for (k in 1:no_treat){

index[k] = Gittins_mat[s[k], f[k]]

}

z[j,1] = covariatedata[j] 

for (v in 1:no_cat){

if (z[j, (t+1)] == v){
          
          max_index = max(index[seq(from = v, to = no_treat - no_cat + v, by = no_cat)])
          
          S = which(index[seq(from = v, to = no_treat - no_cat + v, by = no_cat)] == max_index)
          
          if (length(S) > 1){ 
          
          posi = floor(runif(1, 1, length(S) + 1))
          
          kmax = (S[posi] - 1)*no_cat + v
          
          }else{
          
          kmax = (S-1)*no_cat + v }
          
          } 
        
        }
        
      selected[j, t+1] = kmax
      
      prob_suc_kmax = s[kmax]/(s[kmax] + f[kmax]) 
      
      Pos = as.numeric(runif(1) <= prob_suc_kmax)
      
      if (Pos == 1){
      
      s[kmax] = s[kmax] + 1
      
      }else{
      
      
      s[kmax] = f[kmax] + 1 }
      
      z[j, t+2] = round(runif(1, min = 1, max = no_cat))
      
      selected[j, block+1] = z[j,1]
    
      } 
    
    }
  
  pi_comb_arms = matrix(0L, noRuns, no_treat + 1)
  
  for (j in 1:noRuns){
  
  for (k in 1:no_treat){
  
  pi_comb_arms[j,k] = sum(selected[j, 1:block] == k)
  
  }
  
  pi_comb_arms[j, no_treat + 1] = selected[j, block + 1]
  
  }
  
  id = list()
  
  for (l in 1:no_cat){
  
  id[[l]] = which(pi_comb_arms[ ,no_treat + 1] == l)
  
  }
  
  Totprobmat = matrix(NA, nrow = no_cat, ncol = no_treat)
  
  for (l in 1:no_cat){
  
  if(length(id[[l]]) > 1){
  
  Totprobmat[l,] = colMeans(pi_comb_arms[id[[l]], 1:no_treat])*mean(covariatedata == l)
  
  }else if(length(id[[l]]) == 1){
  
  Totprobmat[l,] = (pi_comb_arms[id[[l]], 1:no_treat])*mean(covariatedata == l)
  
  }else if(length(id[[l]]) == 0){
  
  Totprobmat[l,] = rep(0, no_treat) }
  
  }
  
  Totprob = colSums(Totprobmat)
  
  allocation_probabilities_list = list()
  
  for (l in 1:no_cat){
  
  if(block == 1){
  
  allocation_probabilities_list[[l]] = (Totprob[seq(from = l, to = no_treat - no_cat + l, by = no_cat)])/(mean(z[,1] == l))
  
  }else{
  
  allocation_probabilities_list[[l]] = (Totprob[seq(from = l, to = no_treat - no_cat + l, by = no_cat)])/(mean(rowMeans(z[,1:block] == l)*block))}
 
 }
 
 allocation_probabilities = matrix(unlist(allocation_probabilities_list), nrow = no_cat, ncol = no_treat/no_cat, byrow = TRUE)
 
 allocation_probabilities

}

# Code to run FLGI Simulations

FLGI = function(k, sample_size, block, no_treat, prob_vec, prior_mat){
  

set.seed(k)

if (floor((sample_size/block))*block != sample_size){

stop("Block Size must divide sample size")

}

if (length(prob_vec) != no_treat){

stop("Number of treatments inconsistent")

}

allocation_prob = NULL

a = matrix(0L, 1, block)

treatment = c()

outcome = c()

no_s = matrix(0L, no_treat, 1) 

no_f = matrix(0L, no_treat, 1)

results = matrix(0L, nrow = no_treat, ncol = 2)


for (t in 0:((sample_size/block) - 1)){

alloc_p = matrix(0L,1 ,no_treat)

I = matrix(c(no_s[,1] + 1, no_f[,1] + 1), byrow = FALSE, ncol = 2)

alloc_p = Input_alloc(I, block, 100, 1, no_treat)

allocation_prob = rbind(allocation_prob, alloc_p)

alloc_p = matrix((c(cumsum(c(0, alloc_p[])))), byrow = TRUE, nrow = 1)


Pob = c()

Pos = c()

for (p in 1:block){

Pob[p] = runif(1)

for (k in 1:no_treat){

if ((Pob[p] > alloc_p[1,k]) & (Pob[p] <= alloc_p[1,k+1])){

a[p] = k

treatment[p + (t*block)] = k - 1

w = runif(1)

Pos[p] = as.numeric(w <= prob_vec[k])

if (Pos[p] == 1){

no_s[k,1] = no_s[k,1] + 1 

outcome[p + (t*block)] = 1

}else{

no_f[k,1] = no_f[k,1] + 1

outcome[p + (t*block)] = 0 } 
      }

    } 

  }

}


results[,1] = no_s[,1]

results[,2] = no_f[,1]

outcome = c(results[,1], results[,2])


p.values = rep(NA, no_treat-1)

for (i in 2:no_treat) {

if(nrow(results) == 2){

test.results = results

}else{

test.results = rbind(results[1,], results[i,])}

test.res = fisher.test(test.results, alternative = "less")

p.values[i-1] = test.res$p.value

}

output = list("success" = no_s, "failure" =  no_f, "results" = results, "p.values" = p.values)

return(c(outcome, p.values))

}


# Code to run CARA Equal Randomisation Simulations

CARA_ER = function(k,sample_size, block, no_treat, prob_mat, q){


set.seed(k)

no_s = matrix(0L, ncol = no_treat, nrow = 2)

no_f = matrix(0L, ncol = no_treat, nrow = 2)

cov_pos = vector()

treatment = vector()

outcome = vector()

results_0 = matrix(0L, nrow = no_treat, ncol = 2)

results_1 = matrix(0L, nrow = no_treat, ncol = 2)


for (t in 0:((sample_size/block) - 1)) {

for(p in 1:block){

cov_pos[p + (t*block)] = sample(c(1,0), size = 1, prob = c(q, 1-q), replace = TRUE)

if(cov_pos[p + (t*block)] == 0){

treatment[p + (t*block)] = sample(no_treat, 1, prob = rep(1/no_treat, no_treat), replace = TRUE) # equal probability of being assigned

}else{

treatment[p + (t*block)] = sample(no_treat, 1, prob = rep(1/no_treat, no_treat), replace = TRUE) }

}


for(p in 1:block){

for (k in 1:no_treat){

if(treatment[p + (t*block)] == k ){

outcome[p + (t*block)] = sample(c(1,0), size = 1, replace = TRUE, prob = c(prob_mat[cov_pos[p + (t*block)] + 1, k], 1 - prob_mat[cov_pos[p + (t*block)] + 1, k]))

  } 

}

if(outcome[p + (t*block)] == 1){

no_s[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] = no_s[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] + 1

}else{

no_f[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] = no_f[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] + 1} 
  }

}


results_0[,1] = no_s[1,] # covariate negative

results_0[,2] = no_f[1,] # covariate negative

results_1[,1] = no_s[2,] # covariate positive

results_1[,2] = no_f[2,] # covariate positive

Negative_outcome = c(results_0[1,1], results_0[2,1], results_0[1,2], results_0[2,2])

Positive_outcome = c(results_1[1,1], results_1[2,1], results_1[1,2], results_1[2,2])


trt = vector()

for(t in 1:sample_size){

trt[t] = treatment[t] - 1

} # -1 is to make treatments between 0 and 1

Table = data.frame(outcome, trt, cov_pos)

Model = glm(outcome ~ trt + cov_pos + trt:cov_pos, family = binomial, data = Table)

coef_vec = coef(Model) # extracting model coefficients

pval_vec = coef(Model) #sets up p value vector with NA in right place

pval_vec[!is.na(pval_vec)] = coef(summary(Model))[,4] #replaces non-NA elements with actual p values

#function to merge the two vectors

merge = function(v1,v2){

out = c()

out[seq(from=1,by=2,length.out = length(v1))] = v1

out[seq(from=2,by=2,length.out = length(v2))] = v2

return(out)

}

Model_outcome =  merge(coef_vec,pval_vec)

All = c(Negative_outcome, Positive_outcome, Model_outcome)

return(All)

}


# Code to run CARA Thompson Sampling Simulations

CARA_TS = function(k, sample_size, block, no_treat, prob_mat, prior_mat0, prior_mat1, q){
  
  
  set.seed(k)
  
  no_s = matrix(0L, ncol = no_treat, nrow = 2)
  
  no_f = matrix(0L, ncol = no_treat, nrow = 2)
  
  cov_pos = vector()
  
  treatment = vector()
  
  outcome = vector()
  
  results_0 = matrix(0L, nrow = no_treat, ncol = 2)
  
  results_1 = matrix(0L, nrow = no_treat, ncol = 2)
  
  
  for (t in 0:((sample_size/block) - 1)) {
  
  #Calculating the allocation probabilities for each treatment/covariate combination at the start of each block
  
  theta_01 = rbeta(100, no_s[1,2] + prior_mat0[2,1], no_f[1,2] + prior_mat0[2,2])
  
  theta_00 = rbeta(100, no_s[1,1] + prior_mat0[1,1], no_f[1,1] + prior_mat0[1,2])
  
  theta_neg = cbind(theta_01, theta_00)
  
  allocp_neg = sum(theta_neg[,1] > theta_neg[,2])/(nrow(theta_neg))
  
  
  theta_11 = rbeta(100, no_s[2,2] + prior_mat1[2,1], no_f[2,2] + prior_mat1[2,2])
  
  theta_10 = rbeta(100, no_s[2,1] + prior_mat1[1,1], no_f[2,1] + prior_mat1[1,2])
  
  theta_pos = cbind(theta_11, theta_10)
  
  allocp_pos = sum(theta_pos[,1] > theta_pos[,2])/(nrow(theta_pos))
  
  
  for(p in 1:block){
  
  cov_pos[p + (t*block)] = sample(c(1,0), size = 1, prob = c(q, 1-q), replace = TRUE)
  
  if(cov_pos[p + (t*block)] == 0)
  
  treatment[p + (t*block)] = sample(no_treat, size = 1, replace = TRUE, prob = c(1-allocp_neg, allocp_neg))
  
  else{
  
  treatment[p + (t*block)] = sample(no_treat, size = 1, replace = TRUE, prob = c(1-allocp_pos, allocp_pos)) }
  
  for (k in 1:no_treat){
  
  if(treatment[p + (t*block)] == k ){
  
  outcome[p + (t*block)] = sample(c(1,0), size = 1, replace = TRUE, prob = c(prob_mat[cov_pos[p + (t*block)] + 1, k], 1 - prob_mat[cov_pos[p + (t*block)] + 1, k]))}
  
  }
  
  if(outcome[p + (t*block)] == 1){
  
  no_s[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] = no_s[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] + 1
  
  }else{
  
  no_f[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] = no_f[cov_pos[p + (t*block)] + 1, treatment[p + (t*block)]] + 1} 
    }
  
  }
  
  
  results_0[,1] = no_s[1,] # covariate negative
  
  results_0[,2] = no_f[1,] # covariate negative
  
  results_1[,1] = no_s[2,] # covariate positive
  
  results_1[,2] = no_f[2,] # covariate positive
  
  Negative_outcome = c(results_0[1,1], results_0[2,1], results_0[1,2], results_0[2,2])
  
  Positive_outcome = c(results_1[1,1], results_1[2,1], results_1[1,2], results_1[2,2])
  
  
  trt = vector()
  
  for(t in 1:sample_size){
  
  trt[t] = treatment[t] - 1
  
  } # -1 is to make treatments between 0 and 1
  
  
  Table = data.frame(outcome, trt, cov_pos)
  
  Model = glm(outcome ~ trt + cov_pos + trt:cov_pos, family = binomial, data = Table)
  
  
  coef_vec = coef(Model) # extracting model coefficients
  
  pval_vec = coef(Model) #sets up p value vector with NA in right place
  
  pval_vec[!is.na(pval_vec)] = coef(summary(Model))[,4] #replaces non-NA elements with actual p values
  
  
  #function to merge the two vectors
  
  merge = function(v1,v2){
  
  out = c()
  
  out[seq(from=1,by=2,length.out = length(v1))] = v1
  
  out[seq(from=2,by=2,length.out = length(v2))] = v2
  
  return(out)
  
  }
  
  
  Model_outcome =  merge(coef_vec,pval_vec)
  
  All = c(Negative_outcome, Positive_outcome, Model_outcome)
  
  return(All)
}


# Code to run CARA Forward-Looking Gittins Index Simulations

# Gittins_mat and Input_alloc are the same as in the non-adjusted FLGI code 

CARA_FLGI = function(k, sample_size, block, no_treat, prob_mat, prior_mat0, prior_mat1, q){


set.seed(k)

if (floor((sample_size/block))*block != sample_size){

stop("Block Size must divide sample size")

}


cov_pos = c()

allocation_prob = NULL

no_s = matrix(0L, ncol = no_treat, nrow = 2)

no_f = matrix(0L, ncol = no_treat, nrow = 2)

a = matrix(0L, 1, block)

treatment = c()

outcome = c()

results_0 = matrix(0L, nrow = no_treat, ncol = 2)

results_1 = matrix(0L, nrow = no_treat, ncol = 2)


for (t in 0:((sample_size/block) - 1)){

for (p in 1:block){

cov_pos[p + (t*block)] = sample(c(1,0), size = 1, prob = c(q, 1-q), replace = TRUE)

}


alloc_p = matrix(0L,1 ,no_treat)

for(i in 1:no_treat){

if(cov_pos[p + (t*block)] == 0){

I0 = matrix(c(no_s[1,] + prior_mat0[i,1] + 1, no_f[1,] + prior_mat0[i,2] + 1), byrow = FALSE, ncol = 2)

alloc_p = Input_alloc(I0, block, 100, 1, no_treat)

}else{

I1 = matrix(c(no_s[2,] + prior_mat1[i,1] + 1, no_f[2,] + prior_mat1[i,2] + 1), byrow = FALSE, ncol = 2)

alloc_p = Input_alloc(I1, block, 100, 1, no_treat)

  }

}

allocation_prob = rbind(allocation_prob, alloc_p)

alloc_p = matrix((c(cumsum(c(0, alloc_p[])))), byrow = TRUE, nrow = 1)


Pob = c()

Pos = c()


for (p in 1:block){


Pob[p] = runif(1)

for (k in 1:no_treat){

if ((Pob[p] > alloc_p[1,k]) & (Pob[p] <= alloc_p[1, k+1])){

a[p] = k

treatment[p + (t*block)] = k - 1

w = runif(1)

Pos[p] = as.numeric(w <= prob_mat[cov_pos[p] + 1, k])

if (Pos[p] == 1){

no_s[ cov_pos[p + (t*block)] + 1,k] = no_s[ cov_pos[p + (t*block)] + 1,k] + 1

outcome[p + (t*block)] = 1

}else{

no_f[ cov_pos[p + (t*block)] + 1,k] = no_f[ cov_pos[p + (t*block)] + 1,k] + 1

outcome[p + (t*block)] = 0 } 

      } 

    } 

  } 

}


results_0[,1] = no_s[1,] # covariate negative

results_0[,2] = no_f[1,] # covariate negative

results_1[,1] = no_s[2,] # covariate positive

results_1[,2] = no_f[2,] # covariate positive

Negative_outcome = c(results_0[1,1], results_0[2,1], results_0[1,2], results_0[2,2])

Positive_outcome = c(results_1[1,1], results_1[2,1], results_1[1,2], results_1[2,2])


Table = data.frame(outcome, treatment, cov_pos)

Model = glm(outcome ~ treatment + cov_pos + treatment:cov_pos, family = binomial, data = Table)


coef_vec = coef(Model) # extracting model coefficients

pval_vec = coef(Model) #sets up p value vector with NA in right place

pval_vec[!is.na(pval_vec)] = coef(summary(Model))[,4] #replaces non-NA elements with actual p values


#function to merge the two vectors

merge = function(v1,v2){

out = c()

out[seq(from=1,by=2,length.out = length(v1))] = v1

out[seq(from=2,by=2,length.out = length(v2))] = v2

return(out)

}


Model_outcome =  merge(coef_vec,pval_vec)

All = c(Negative_outcome, Positive_outcome, Model_outcome)

return(All)

}
