# Dissertation
All code to run the simulations in my Masters Dissertation

Code to run Equal Randomisation simulations

Equal_random = function(k, sample_size, block, no_treat, prob_vec){
  
  set.seed(k)
  
  treatment = vector()
  
  outcome = vector()
  
  no_s = matrix(0L, no_treat, 1)
  
  no_f = matrix(0L, no_treat, 1)
  
  results = matrix(0L, nrow = no_treat, ncol = 2)
  
  
  for (t in 0:((sample_size/block) - 1)) {
    
    treatment = sample(no_treat, sample_size, prob = rep(1/no_treat, no_treat), replace = TRUE) # equal probability of being assigned
    
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
  
  p_val = p.adjust(p.values, method = "bonferroni")
  
  output = list("success" = no_s, "failure" =  no_f, "results" = results, "p.values" = p.values)
  
  return(c(outcome, p_val))

}

