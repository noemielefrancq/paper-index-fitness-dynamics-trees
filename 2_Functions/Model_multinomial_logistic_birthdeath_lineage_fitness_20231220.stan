//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//
functions{
  vector get_alpha(vector alpha_true, array[] int parents, int K, int G, vector t_start, vector t){
    vector[K] alpha=rep_vector(-1,K);   // Offset ancestral lineage for each structural group (-1 is group is not ancestral)
    int count_alpha;
    count_alpha=1;
    for (k in 1:K){
      if (parents[k]==0 && t_start[k] <= t[1]){ //if ancestral group present from the begenning
        alpha[k] = alpha_true[count_alpha];
        if(alpha[k] == 0) alpha[k] = 1E-30; // To avoid model crash
        count_alpha = count_alpha+1;
     }
   }
   return(alpha);
  }
  
  vector get_alpha_GA(vector alpha_true_GA, array[] int parents, int K, int G, vector t_start, vector t){
    vector[K] alpha_GA=rep_vector(-1,K);   // Offset ancestral lineage for each structural group (-1 is group is not ancestral)
    int count_alpha;
    count_alpha=1;
    for (k in 1:K){
      if (parents[k]==0 && t_start[k] > t[1]){ //if ancestral group present from the begenning
        alpha_GA[k] = alpha_true_GA[count_alpha];
        if(alpha_GA[k] == 0) alpha_GA[k] = 1E-30; // To avoid model crash
        count_alpha = count_alpha+1;
     }
   }
   return(alpha_GA);
  }
  
  int N_unique(array[] int x){
    array[num_elements(x)] int x_sorted = sort_asc(x);
    array[num_elements(x)] int tmp = rep_array(1,num_elements(x));
    array[num_elements(x)] int tmp_idx = rep_array(1,num_elements(x));
    for (i in 2:num_elements(x_sorted)) {
      if(x_sorted[i-1] == x_sorted[i]){
        tmp[i] = tmp[i-1] + 1;
      }
    }
    tmp_idx = sort_indices_asc(tmp);
    int num = 0;
    for(i in 1:num_elements(tmp_idx)){
      if(tmp_idx[i] == 1) num = num + 1;
    }
    return(num);
  }
  
  array[] int unique(array[] int x){
    array[num_elements(x)] int x_sorted = sort_asc(x);
    array[num_elements(x)] int tmp = rep_array(1,num_elements(x));
    array[num_elements(x)] int tmp_idx = rep_array(1,num_elements(x));
    for (i in 2:num_elements(x_sorted)) {
      if(x_sorted[i-1] == x_sorted[i]){
        tmp[i] = tmp[i-1] + 1;
      }
    }
    tmp_idx = sort_indices_asc(tmp);
    int num = 0;
    for(i in 1:num_elements(tmp_idx)){
      if(tmp_idx[i] == 1) num = num + 1;
    }
    array[num] int res = x_sorted[tmp_idx[1:num]];
    return(res);
  }
  
  vector get_gamma(vector gamma_true, array[] int parents, int K){
    vector[K] gamma=rep_vector(0,K);   //share of the parent group in the new group (0 if group is ancestral)
    int count_gamma;
    count_gamma=1;
    for (k in 1:K){
      if (parents[k]!=0){ //if not ancestral group
        gamma[k]=gamma_true[count_gamma];
        count_gamma=count_gamma+1;
     }
   }
   return(gamma);
  }
  
  array[] int get_lin_presence(int n, int num_lin, array[,] int lin_presence){
    array[num_lin] int which_lin_presence;
    which_lin_presence = rep_array(0, num_lin);
    int count_lineage;
    count_lineage=1;
    int id_lineage;
    id_lineage=1;
    for (l in 1:dims(lin_presence)[2]){
      if (lin_presence[n,l] == 1){
        which_lin_presence[count_lineage] = id_lineage;
        count_lineage+=1;
      }
      id_lineage+=1;
    }
    return(which_lin_presence);
  }
  
  int num_matches(array[] int x, int a) {
    int n = 0;
    for (i in 1:size(x))
      if (x[i] == a)
        n += 1;
    return n;
  }
  
  int num_below(vector x, real a) {
    int n = 0;
    for (i in 1:size(x))
      if (x[i] <= a)
        n += 1;
    return n;
  }
  
  array[] int which_equal(array[] int x, int a) {
    array[num_matches(x, a)] int match_positions;
    int pos = 1;
    for (i in 1:size(x)) {
      if (x[i] == a) {
        match_positions[pos] = i;
        pos += 1;
      }
    }
    return match_positions;
  }
  
  array[] int which_below(vector x, real a) {
    array[num_below(x, a)] int match_positions;
    int pos = 1;
    for (i in 1:size(x)) {
      if (x[i] <= a) {
        match_positions[pos] = i;
        pos += 1;
      }
    }
    return match_positions;
  }
  
  matrix compute_theta(vector alpha, vector alpha_GA, vector beta, vector gamma, array[] int parents, array[,] int lin_presence, vector t_start, vector t, int K, int N){
    matrix[K, N] theta = rep_matrix(-1e30, K, N);
    int num_lin_init = num_below(t_start, t[1]);
    array[num_lin_init] int which_lin_presence_init;
    which_lin_presence_init = which_below(t_start, t[1]);
    theta[which_lin_presence_init,1] = log(alpha[which_lin_presence_init]);//in the beginning, only ancestral lineages are present
    for (i in which_lin_presence_init){
      if (alpha[i] == -1){//specific case where a new lineage appears on the first time step
        theta[i, 1] = log(gamma[i]) + theta[parents[i],1];
      }
    }
    if(num_lin_init > 1){
      theta[,1] = theta[,1] - theta[K, 1]; // Make sure the reference is at 0
    } // Else: reference is already at 0 
    for (n in 2:N){  //for each year
    
      // Get IDs of the lineages that were present at time t-1
      int num_lin_t1 = num_below(t_start, t[n-1]); 
      array[num_lin_t1] int which_lin_presence_t1;
      which_lin_presence_t1 = which_below(t_start, t[n-1]);

      // Get IDs of the lineages that are present at time t
      int num_lin_t = num_below(t_start, t[n]);
      array[num_lin_t] int which_lin_presence_t;
      which_lin_presence_t = which_below(t_start, t[n]);

      array[K] int is_new_lineage;
      is_new_lineage = rep_array(-10,K);
      for (k in which_lin_presence_t){
        if (k != K){
          if (num_matches(which_lin_presence_t1, k) == 1){ //if the lineage was here before
            theta[k,n] = theta[k,(n-1)]+beta[k]*(t[n]-t[n-1]);
          }
          else{ //if new lineage
            is_new_lineage[k] = parents[k]; // Record the parent of each new lineage
          }
        }else{
          theta[k,n] = theta[k,(n-1)]; // reference
        }
      }
      for (i in 1:K){ // We deal with parent sepratelty, and make the new lienages appear 
        int n_parents = num_matches(is_new_lineage,  K-i+1); // Count number of parents that produce a new lieange at this time step
        if(n_parents > 0){
          // Get indexes of the new lineages
          array[n_parents] int idx;
          idx = which_equal(is_new_lineage, K-i+1);
          idx = idx[sort_indices_asc(t_start[idx])]; // sort idx by increaseing t_starts
          
          // Store the parent initial theta (at time t-1, or t_start if the parent also appeared wihtin this time step)
          real tmp;
          if(num_matches(which_lin_presence_t1, K-i+1) == 1){
            if(K-i+1 != K){
              tmp = theta[K-i+1,(n-1)]+beta[K-i+1]*(t_start[idx[1]]-t[n-1]); // if parent was there at time n-1
            }else{
              tmp = theta[K-i+1,(n-1)]; // parent is the reference group
            }
          }else{
            tmp = theta[K-i+1,n] - beta[K-i+1]*(t_start[idx[1]]-t_start[K-i+1]);
          }
          
          // Rather complex loop to compute the thetas of the offspring, and update parent's theta
          for (j in 1:num_elements(idx)){
            theta[idx[j], n] = log(gamma[idx[j]]) + tmp + beta[idx[j]]*(t[n] - t_start[idx[j]]);
            if(num_elements(idx) == j){
              if(K-i+1 != K){
                theta[K-i+1, n] = log(1-gamma[idx[j]]) + tmp + beta[K-i+1]*(t[n] - t_start[idx[j]]);
              }else{
                theta[K-i+1, n] = tmp; // parent is the reference group, by definition reference group's theta remains constant
              }
            }else{
              if(K-i+1 != K){
                tmp = log(1-gamma[idx[j]]) + tmp + beta[K-i+1]*(t_start[idx[j+1]] - t_start[idx[j]]);
              }
            }
          }
        }
      }
      // Potential case of 'new lienage but no parent'
      int n_parents = num_matches(is_new_lineage,  0); // Count number of parents that produce a new lieange at this time step
      if(n_parents > 0){
        // Get indexes of the new lineages
        array[n_parents] int idx;
        idx = which_equal(is_new_lineage, 0);
        
        real tmp = log(sum(exp(theta[, n])));
        for (j in 1:num_elements(idx)){
          theta[idx[j], n] = log(alpha_GA[idx[j]]) + tmp;
        }
      }
    }
    return theta;
  }
  
  real multinomial_zeros_lpmf(array[] int y, vector p){
    // y contains the counts per category
    // p is a simplex: probabilities of observing heach category
    real res;
    
    // Usal multinomial
    res = multinomial_lpmf(y | p);
    
    // Count the number of categories with 0 counts: if any, assess the likelihood of observing 0 counts given the predicted proportion of the group in the population
    int n_zeros_in_y = num_matches(y, 0);
    if(n_zeros_in_y > 0){
      array[n_zeros_in_y] int idx;
      idx = which_equal(y, 0);
      for (i in 1:n_zeros_in_y){
        res = res + binomial_lpmf(0 | sum(y), p[idx[i]]); // Likelihoof of a category with 0 count
      }
    }
  
    return(res);
  }
}

data{
  // Data
  int <lower = 0> N;       // Number of time points
  int <lower=0> G;         // Number of structural groups, present from the begenning
  int <lower=0> GA;        // Number of structural groups, that appear at some point
  int <lower = 0> K;       // Total number of lineages
  
  array[N,K] int<lower=0> Y;     // Number of sequences of each lineage at each time point
  vector[N] t;             // Time
  array[K] int <lower=0> parents;       //parents of each lineage (0 if no parents)
  vector[K] t_start;       // Origin of linev ages
  array[N,K] int lin_presence; //Presence or absence of each lineage at each time point

  // To produce more detailed plots
  int <lower = 0> N_new;
  vector[N_new] t_new;
}

parameters {
  // Proportions
  simplex[G] alpha_true;   // Starting proportion of ancestral groups 
  vector<lower=0, upper = 1> [GA] alpha_true_GA; // Starting proportion of each new group (ancestral groups that appear after start)
  vector<lower=0, upper = 1> [K-G-GA] gamma_true; // Starting proportion of each new group (not ancestral)
    
  // Fitness
  vector[K-1] beta;        // Growth rate of each group 
}

transformed parameters { 
  // Get alpha for each group
  vector[K] alpha=rep_vector(0, K);
  alpha=get_alpha(alpha_true, parents, K, G, t_start, t);
   
  // Get true alpha GA
  vector[K] alpha_GA=rep_vector(0, K);
  if(GA > 0){
    alpha_GA=get_alpha_GA(alpha_true_GA, parents, K, G, t_start, t);
  }
  
  // Get true gamma
  vector[K] gamma=rep_vector(0, K);
  gamma=get_gamma(gamma_true, parents, K);
  
  // Theta theta
  matrix[K, N] theta = rep_matrix(-1e30, K, N);
  theta=compute_theta(alpha, alpha_GA, beta, gamma, parents, lin_presence, t_start, t, K, N);
}

model {
   beta ~ double_exponential(0, 1);
   alpha_true_GA ~ normal(0, 0.05);
   // gamma_true ~ normal(0, 0.05); // Prior on starting frequencies of new groups: we expect small frequencies
   gamma_true ~ beta(1, 99); // Prior on starting frequencies of new groups: we expect small frequencies

   for (n in 1:N){
     vector[K] tmp = softmax(theta[, n]);
     for(i in 1:K){
       if(tmp[i] < 1e-10) tmp[i] = 1e-10; // Set minimum frequency to 0.00001 (to avoid log(prop) = -inf)
     }
     tmp = tmp/sum(tmp); // Renormalize
     
     // Likelihood
     Y[n,] ~ multinomial_zeros(tmp);
   }
}

generated quantities {
  // Compute likelihood of each data point (useful for model comparison)
  vector[N] log_lik;
  for (n in 1:N) {
     vector[K] tmp = softmax(theta[, n]);
     for(i in 1:K){
       if(tmp[i] < 1e-5) tmp[i] = 1e-5;
     }
     tmp = tmp/sum(tmp);
     log_lik[n] = multinomial_zeros_lpmf(Y[n,] | tmp);
  }
  
  // Compute frequencies at more time points (to have better plots)
  matrix[K, N_new] theta_new = rep_matrix(-1e30, K, N_new);
  theta_new = compute_theta(alpha, alpha_GA, beta, gamma, parents, lin_presence, t_start, t_new, K, N_new);
}



