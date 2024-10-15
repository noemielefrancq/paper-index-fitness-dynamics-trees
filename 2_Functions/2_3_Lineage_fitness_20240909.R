## Codes for fitness estimation

## Tool functions
mean.and.ci <-function(v){ return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}
read.chains.from.table = function(table){
  if(typeof(table) != 'double') {
    print('Changing type of table to matrix')
    table = as.matrix(table)
  }
  Chains = list()
  col_variables = sapply(colnames(table), function(x)str_split(x, pattern = "[.]")[[1]][1])
  variable_names = unique(col_variables)
  nchains = nrow(table)
  for(i in 1:length(variable_names)){
    a = match(col_variables, variable_names[i])
    if(length(which(is.na(a) == F)) == 1) {
      Chains[[i]] = table[,match(variable_names[i], col_variables)]
    }
    else{
      a = match(col_variables, variable_names[i])
      tmp = colnames(table)[which(is.na(a) == F)]
      ndims = 0
      dims = NULL
      empty = F
      while(empty == F & ndims < 10){
        l = sapply(tmp, function(x)str_split(x, pattern = "[.]")[[1]][1+ndims+1])
        if(all(is.na(l))) empty = T
        if(all(is.na(l)) == F) {
          ndims = ndims + 1
          dims = c(dims, max(as.numeric(l)))
        }
      }
      if(ndims >5) print('Error: this function only supports arrays of <=5 dimensions')
      Chains[[i]] = array(NA, dim = c(nchains, dims))
      a = which(is.na(match(col_variables, variable_names[i]))==F)
      
      if(ndims == 1){
        Chains[[i]] = table[,a]
        colnames(Chains[[i]]) = NULL
      }
      
      if(ndims == 2){
        j = 1
        for(d2 in 1:dims[2]){
          for(d1 in 1:dims[1]){
            Chains[[i]][,d1,d2] = table[,a[j]]
            j = j+1
          }
        }
      }
      
      if(ndims == 3){
        j = 1
        for(d3 in 1:dims[3]){
          for(d2 in 1:dims[2]){
            for(d1 in 1:dims[1]){
              Chains[[i]][,d1,d2,d3] = table[,a[j]]
              j = j+1
            }
          }
        }
      }
      
      if(ndims == 4){
        j = 1
        for(d4 in 1:dims[4]){
          for(d3 in 1:dims[3]){
            for(d2 in 1:dims[2]){
              for(d1 in 1:dims[1]){
                Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                j = j+1
              }
            }
          }
        }
      }
      
      if(ndims == 5){
        j = 1
        for(d5 in 1:dims[5]){
          for(d4 in 1:dims[4]){
            for(d3 in 1:dims[3]){
              for(d2 in 1:dims[2]){
                for(d1 in 1:dims[1]){
                  Chains[[i]][,d1,d2,d3,d4] = table[,a[j]]
                  j = j+1
                }
              }
            }
          }
        }
      }
    }
    names(Chains)[i] = variable_names[i]
  }
  return(Chains)
}

## Functions to plot fits
plot_fit_data = function(data, Chains, colour_lineage){
  plot(NULL, bty = "n", ylim = c(0,1), xlim = c(min(data$t), max(data$t)), xlab = "Time (years)", ylab = "Proportion")
  pred_freq_chains = array(NA, dim = c(length(Chains$lp__), data$K, data$N))
  for(i in 1:data$N){
    pred_freq_chains[,,i] = t(apply(Chains$theta[,,i], MAR = 1, softmax))
  }
  data_m = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimin = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimax = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  tmp = lapply(1:data$N, function(x)DescTools::MultinomCI(data$Y[x,]))
  for(i in 1:data$N){
    data_m[i,] = tmp[[i]][,1]
    data_cimin[i,] = tmp[[i]][,2]
    data_cimax[i,] = tmp[[i]][,3]
  }
  for(i in 1:data$K){
    pred_freq = apply(pred_freq_chains[,i,], MARGIN = 2, function(x)mean.and.ci(x))
    ## fit
    lines(data$t, pred_freq[1,], lwd = 2, col = colour_lineage[i])
    polygon(x = c(data$t, rev(data$t)),
            y = c(pred_freq[2,], rev(pred_freq[3,])), border = F,
            col = adjustcolor(colour_lineage[i], alpha.f = 0.5))
    ## data, Multinomial CI
    d_m = data_m[,i]
    d_cimin = data_cimin[,i]
    d_cimax = data_cimax[,i]
    points(data$t, d_m, col = colour_lineage[i], pch = 16)
    arrows(data$t, d_cimin, data$t, d_cimax, length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor(colour_lineage[i], alpha.f = 0.8))
  }
}
plot_fit_data_new = function(data, Chains, colour_lineage, xmin, xmax){
  plot(NULL, bty = 'n', ylim = c(0,1), xlim = c(xmin, xmax), yaxt = 'n',
       col = colour_lineage[1], xlab = 'Time (years)', ylab = 'Proportion')
  axis(2, las = 2)
  pred_freq_chains = array(NA, dim = c(length(Chains$lp__), data$K, data$N_new))
  for(i in 1:data$N_new){
    pred_freq_chains[,,i] = t(apply(Chains$theta_new[,,i], MAR = 1, softmax))
  }
  data_m = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimin = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimax = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  tmp = lapply(1:data$N, function(x)DescTools::MultinomCI(data$Y[x,]))
  for(i in 1:data$N){
    data_m[i,] = tmp[[i]][,1]
    data_cimin[i,] = tmp[[i]][,2]
    data_cimax[i,] = tmp[[i]][,3]
  }
  for(i in 1:data$K){
    pred_freq = apply(pred_freq_chains[,i,], MARGIN = 2, function(x)mean.and.ci(x))
    ## fit
    lines(data$t_new, pred_freq[1,], lwd = 2, col = colour_lineage[i])
    polygon(x = c(data$t_new, rev(data$t_new)),
            y = c(pred_freq[2,], rev(pred_freq[3,])), border = F,
            col = adjustcolor(colour_lineage[i], alpha.f = 0.5))
    ## data, Multinomial CI
    d_m = data_m[,i]
    d_cimin = data_cimin[,i]
    d_cimax = data_cimax[,i]
    ## Plot
    points(data$t, d_m, col = colour_lineage[i], pch = 16)
    arrows(data$t, d_cimin, data$t, d_cimax, length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor(colour_lineage[i], alpha.f = 0.8))
  }
}
plot_fit_data_per_group = function(data, Chains, colour_lineage){
  pred_freq_chains = array(NA, dim = c(length(Chains$lp__), data$K, data$N))
  for(i in 1:data$N){
    pred_freq_chains[,,i] = t(apply(Chains$theta[,,i], MAR = 1, softmax))
  }
  data_m = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimin = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimax = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  tmp = lapply(1:data$N, function(x)DescTools::MultinomCI(data$Y[x,]))
  for(i in 1:data$N){
    data_m[i,] = tmp[[i]][,1]
    data_cimin[i,] = tmp[[i]][,2]
    data_cimax[i,] = tmp[[i]][,3]
  }
  for(i in 1:data$K){
    plot(NULL, bty = "n", ylim = c(0,1), xlim = c(min(data$t), max(data$t)), xlab = "Time (years)", ylab = "Proportion",
         main = i)
    pred_freq = apply(pred_freq_chains[,i,], MARGIN = 2, function(x)mean.and.ci(x))
    ## fit
    lines(data$t, pred_freq[1,], lwd = 2, col = colour_lineage[i])
    polygon(x = c(data$t, rev(data$t)),
            y = c(pred_freq[2,], rev(pred_freq[3,])), border = F,
            col = adjustcolor(colour_lineage[i], alpha.f = 0.5))
    ## data, Multinomial CI
    d_m = data_m[,i]
    d_cimin = data_cimin[,i]
    d_cimax = data_cimax[,i]
    ## Plot
    points(data$t, d_m, col = colour_lineage[i], pch = 16)
    arrows(data$t, d_cimin, data$t, d_cimax, length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor(colour_lineage[i], alpha.f = 0.8))
  }
}
plot_fit_data_selected = function(data, Chains, colour_lineage, selected){
  plot(NULL, bty = "n", ylim = c(0,1), xlim = c(min(data$t), max(data$t)), xlab = "Time (years)", ylab = "Proportion", yaxt = 'n', xaxt = 'n')
  axis(2, las = 2, lwd = 0.5, lwd.ticks = 0.5, tck=-0.01, mgp = c(0.8, 0.2, 0), at = seq(0,1,0.5), labels = seq(0,1,0.5))
  axis(1, lwd = 0.5, lwd.ticks = 0.5, tck=-0.01, mgp = c(0.8, -0.2, 0))
  pred_freq_chains = array(NA, dim = c(length(Chains$lp__), data$K, data$N_new))
  for(i in 1:data$N_new){
    pred_freq_chains[,,i] = t(apply(Chains$theta_new[,,i], MAR = 1, softmax))
  }
  data_m = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimin = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  data_cimax = matrix(NA, nrow = nrow(data$Y), ncol = ncol(data$Y))
  tmp = lapply(1:data$N, function(x)DescTools::MultinomCI(data$Y[x,]))
  for(i in 1:data$N){
    data_m[i,] = tmp[[i]][,1]
    data_cimin[i,] = tmp[[i]][,2]
    data_cimax[i,] = tmp[[i]][,3]
  }
  for(i in selected){
    pred_freq = apply(pred_freq_chains[,i,], MARGIN = 2, function(x)mean.and.ci(x))
    
    ## Plot only after lineage birth
    t_birth = data$t_start[i]
    
    ## fit
    lines(data$t_new[which(data$t_new >= t_birth)], pred_freq[1,which(data$t_new >= t_birth)], lwd = 2, col = colour_lineage[i])
    polygon(x = c(data$t_new[which(data$t_new >= t_birth)], rev(data$t_new[which(data$t_new >= t_birth)])),
            y = c(pred_freq[2,which(data$t_new >= t_birth)], rev(pred_freq[3,which(data$t_new >= t_birth)])), border = F,
            col = adjustcolor(colour_lineage[i], alpha.f = 0.5))
    ## data, Multinomial CI
    d_m = data_m[which(data$t >= t_birth),i]
    d_cimin = data_cimin[which(data$t >= t_birth),i]
    d_cimax = data_cimax[which(data$t >= t_birth),i]
    ## Plot
    points(data$t[which(data$t >= t_birth)], d_m, col = colour_lineage[i], pch = 16, cex = 0.5)
    arrows(data$t[which(data$t >= t_birth)], d_cimin, data$t[which(data$t >= t_birth)], d_cimax, length=0, angle=0, code=3, lwd = 0.8,
           col = adjustcolor(colour_lineage[i], alpha.f = 0.8))
  }
}
plot_fitness_values = function(data, Chains, colour_lineage, gentime){
  betas = apply(exp(Chains$beta*gentime), MARGIN = 2, function(x)mean.and.ci(x))
  plot(1:(length(betas[1,])+1), c(betas[1,], 1), log = 'y', 
       col = colour_lineage, pch = 16, bty = 'n', ylim = c(0.99,1.3),
       ylab = 'Relative fitness per generation', xaxt = 'n', yaxt = 'n', xlab = '')
  arrows(1:(length(betas[1,])+1), c(betas[3,], 1), 
         1:(length(betas[1,])+1), c(betas[2,], 1), length=0, angle=0, code=3, lwd = 0.8,
         col = colour_lineage)
  abline(h=1, lty = 2)
  axis(2, las = 2)
  axis(side = 1, labels = FALSE, at = 1:(length(betas[1,])+1))
}

## Plot raw fitness estimates
plot_estimated_fitness_ref_ancestral = function(data, Chains, colour_lineage, gentime){
  betas = apply((Chains$beta*gentime), MARGIN = 2, function(x)mean.and.ci(x))
  plot(1:(length(betas[1,])+1), c(betas[1,], 0), cex = 0.5,
       col = colour_lineage, pch = 16, bty = 'n', ylim = c(0,max(betas)),
       ylab = 'Relative fitness', xaxt = 'n', yaxt = 'n', xlab = 'Groups')
  abline(h=0, lty = 2)
  axis(2, las = 2)
  arrows(1:length(betas[1,]), betas[2,], 1:length(betas[1,]), betas[3,], length=0, angle=0, code=3, lwd = 1.5,
         col = adjustcolor(colour_lineage, alpha.f = 0.8))
  axis(side = 1, labels = T, at = 1:(length(betas[1,])+1))
  arrows(1:length(betas[1,]), betas[2,], 1:length(betas[1,]), betas[3,], length=0, angle=0, code=3, lwd = 0.8,
         col = adjustcolor(colour_lineage, alpha.f = 0.8))
}

## Plot observed vs predicted proportions
plot_observed_vs_predicted = function(data, Chains, colour_lineage){
  plot(NULL, bty = 'n', xlim = c(0,1), ylim = c(0,1),
       xlab = 'Observed', ylab = 'Predicted', yaxt = 'n')
  axis(2, las =2)
  abline(b=1, a=0, lty = 2)
  pred_freq_chains = array(NA, dim = c(length(Chains$lp__), data$K, data$N))
  for(i in 1:data$N){
    pred_freq_chains[,,i] = t(apply(Chains$theta[,,i], MAR = 1, softmax))
  }
  for(i in 1:data$K){
    pred_freq = apply((pred_freq_chains[,i,]), MARGIN = 2, function(x)mean.and.ci(x))

    # data
    tmp = binom.confint(x = data$Y[,i], n = rowSums(data$Y), method = c("bayes"), type="central")
    d_m = tmp$mean
    d_cimin = tmp$lower
    d_cimax = tmp$upper
    points(d_m, pred_freq[1,], col = colour_lineage[i], pch = 16, cex = 0.5)
  }
}

## Compute softmax
softmax <- function(x) {exp(x) / sum(exp(x))} 

## Main fit function
estimate_rel_fitness_groups_with_branches = function(dataset_with_nodes, tree, min_year = 1950, window = NULL, N = NULL, 
                                                     model_compiled, iter_warmup = 250, iter_sampling = 500, refresh = 50, seed = 1){
  ## Retrieve branches from tree
  branches = tree$edge
  
  ## Branch times
  branches_times = tree$edge
  branches_times[,1] = branches_times[,2] = NA
  times_tips_nodes = dataset_with_nodes$time
  branches_times[,1] = times_tips_nodes[match(branches[,1], dataset_with_nodes$ID)]
  branches_times[,2] = times_tips_nodes[match(branches[,2], dataset_with_nodes$ID)]
  
  ## Branch groups
  branches_group = tree$edge
  branches_group[,1] = branches_group[,2] = NA
  branches_group[,1] = as.numeric(dataset_with_nodes$groups[branches[,1]])
  branches_group[,2] = as.numeric(dataset_with_nodes$groups[branches[,2]])
  
  ## Group names
  groups = table(factor(dataset_with_nodes$groups))
  groups = names(groups)
  
  ## Set time window, to count number of sequences and nodes within each group
  max_year = max(dataset_with_nodes$time)
  if(is.null(window) == F){
    time_windows = seq(min_year, max(dataset_with_nodes$time), window)
  }else if(is.null(N) == F){
    time_windows=seq(min_year, max_year, length.out=N)
  }
  mid_time = time_windows[-length(time_windows)]+(time_windows[2]-time_windows[1])/2
  
  ## Set data frame to store counts
  count_time = matrix(NA, nrow = length(groups), 
                      ncol = length(mid_time))
  
  ## Compute counts through time
  for(i in 1:length(mid_time)){
    t_min = time_windows[i]
    t_max = time_windows[i+1]
    a = which(branches_times[,1] <= t_min & branches_times[,2] >= t_max) ## branches alive, with no sampled individuals 
    b = which(branches_times[,1] <= t_min & (branches_times[,2] >= t_min & branches_times[,2] < t_max)) ## branches born before interval and died within interval
    c = which((branches_times[,1] >= t_min & branches_times[,1] < t_max) & branches_times[,2] >= t_max) ## branches born in interval and died after interval
    d = which((branches_times[,1] >= t_min & branches_times[,1] < t_max) & (branches_times[,2] >= t_min & branches_times[,2] < t_max) ) ## branches born and died within interval
    
    e = c(a,b,c,d)
    e = e[which(branches_group[e,1] - branches_group[e,2] == 0)] ## Filter, to keep only branches that do not have a switch of group
    
    indiv_time = branches_group[e,1] ## Only consider 1 group per branch
    
    tmp = table(factor(indiv_time, levels = as.numeric(groups))) ## Counts per group
    
    count_time[,i] = tmp # Store results
  }
  
  ## Add information on parents groups
  parents=rep(0, length(as.numeric(groups)))
  for (i in as.numeric(groups)){
    beginning_node=dataset_with_nodes$ID[dataset_with_nodes$groups==i & !is.na(dataset_with_nodes$groups)][1]
    path=nodepath(phy=tree,from=beginning_node, to=tree$Nnode+2)
    group_path=dataset_with_nodes$groups[path]
    unique_path=unique(group_path)
    unique_path=unique_path[!is.na(unique_path)]
    if (length(unique_path)>1) parents[i]=unique_path[2]
  }
  # lineage_pres_abs[parents==0,] = 1 # The ancestral lineages are always here
  # parents[which(lineage_pres_abs[,1] == 1)] = 0 ## Simplify stan run and consider that all the lineages that do exist from the strating point are considered as ancestral in the stan code
  
  ## Compute starting times each group
  t_start = rep(NA, length(groups))
  t_start_upper_bound = rep(NA, length(groups))
  t_start_index = rep(NA, length(groups))
  
  for(j in 1:length(t_start)){
    tmp = which(dataset_with_nodes$groups == groups[j])
    tmp2 = dataset_with_nodes$ID[tmp[which.min(dataset_with_nodes$time[tmp])]]
    if(parents[j] > 0){
      m = branches_times[which(tree$edge[,2] == tmp2),1]
    }else if(parents[j] == 0){
      m = branches_times[which(tree$edge[,2] == tmp2),2]
    }
    if(length(m) == 0) m = min(mid_time) ## Specific case of the root
    t_start[j] = m
    ## Compute index of the starting time (technically not used in the code anymore)
    index = tail(which(mid_time <= t_start[j]), 1)
    if(length(index) == 0) index = 1
    t_start_index[j] = index
  }
  
  ## Compute lineage presence/absence
  lineage_pres_abs =  matrix(0, nrow = length(groups),
                             ncol = length(mid_time))
  for (i in 1:length(groups)){
    lineage_pres_abs[i,t_start_index[i]:dim(lineage_pres_abs)[2]]=1
  }
  parents[which(lineage_pres_abs[,1] == 1)] = 0
  t_start[which(lineage_pres_abs[,1] == 1)] = mid_time[1]-2 ## Make sure the first starting time is before the strat of the time series
  t_start[which.min(t_start)] = mid_time[1]-2
  
  ## Number of ancestral groups that are here from the beginning
  G = length(which(t_start <= mid_time[1] & parents==0))
  
  ## Number of ancestral groups that appear
  GA = length(which(t_start > mid_time[1] & parents==0))
  
  ## Build data object, comprising with all the necessary numbers, vectors and matrices
  data <- list(N=length(mid_time), G=G, GA = GA, K = length(groups), Y = t(count_time),
               parents=parents, 
               t = mid_time, t_start = t_start, t_start_index = t_start_index, 
               t_new = seq(min_year, max_year, length.out = 150), N_new = 150,
               lin_presence = t(lineage_pres_abs))
  
  initial_values = function(){
    return(list('beta' = rnorm(n = length(groups)-1, mean = 0, sd = 0.05), ## Small betas
                'alpha_true' = rmultinom(n=1, size= 100, prob = rep(1/data$G, data$G))[,1]/100, ## Equal-ish starting frequencies of ancestral groups
                'gamma_true' = abs(rnorm(n=length(groups)-data$G, mean=0, sd=0.001)))) ## Small starting frequencies
  }
  
  fit <- model_compiled$sample(data = data, refresh = 50, #seed=24,
                               chains = 3, parallel_chains = 3,
                               iter_warmup = 250, iter_sampling = 250,
                               max_treedepth = 12, adapt_delta = 0.97,
                               init = list(initial_values(), initial_values(), initial_values()))
  
  ## Diagnostic
  fit$cmdstan_diagnose()#check this
  
  ## Extract chains
  t <- list()
  for (f in fit$output_files()) t[[f]] <- as.matrix(data.table::fread(cmd= paste0("grep -v '^#' ", f)))
  t = do.call(rbind, t) ## Combine chains
  t = as.matrix(t)
  Chains = read.chains.from.table(t)
  remove(t)
  
  return(list('fit' = fit,
              'chains' = Chains,
              'data' = data))
}


