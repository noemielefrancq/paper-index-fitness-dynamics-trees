## Codes for finding lineage defining mutations, example on SARS-CoV-2

########################################################################################################################################
## Packages used
########################################################################################################################################
library(stringr)
library(ape)
library(phangorn)
library(BactDating)
library(phytools)
library(coda)
library(thd)
library(vcfR)
library(lubridate)
library(ggplot2)
library(ggtree)
library(extrafont)
library(cowplot)
library(scales)
# font_import()
loadfonts(device="all")
########################################################################################################################################

########################################################################################################################################
## Useful functions
########################################################################################################################################
## Other functions
mean.and.ci <-function(v){ return( c(mean(v), as.numeric(quantile(v,probs = 0.025, na.rm = T)), as.numeric(quantile(v,probs = 0.975, na.rm = T))))}
axisPhylo_NL = function (side = 1, root.time = NULL, backward = TRUE, at_axis = NULL, lab_axis = NULL, ...){
  lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
  type <- lastPP$type
  if (type == "unrooted")
    stop("axisPhylo() not available for unrooted plots; try add.scale.bar()")
  if (type == "radial")
    stop("axisPhylo() not meaningful for this type of plot")
  if (is.null(root.time))
    root.time <- lastPP$root.time
  if (type %in% c("phylogram", "cladogram")) {
    xscale <- if (lastPP$direction %in% c("rightwards", "leftwards"))
      range(lastPP$xx)
    else range(lastPP$yy)
    tmp <- lastPP$direction %in% c("leftwards", "downwards")
    tscale <- c(0, xscale[2] - xscale[1])
    if (xor(backward, tmp))
      tscale <- tscale[2:1]
    if (!is.null(root.time)) {
      tscale <- tscale + root.time
      if (backward)
        tscale <- tscale - xscale[2]
    }
    beta <- diff(xscale)/diff(tscale)
    alpha <- xscale[1] - beta * tscale[1]
    if(is.null(at_axis) == T){
      x <- beta * lab + alpha
      lab <- pretty(tscale)
    }
    # if(is.null(at_axis) != F){
    x <- at_axis
    lab <- lab_axis
    # }
    axis(side = side, at = x, labels = lab, ...)
  }
  else {
    n <- lastPP$Ntip
    xx <- lastPP$xx[1:n]
    yy <- lastPP$yy[1:n]
    r0 <- max(sqrt(xx^2 + yy^2))
    alpha <- sort(setNames(rect2polar(xx, yy)$angle, 1:n))
    angles <- c(diff(alpha), 2 * pi - alpha[n] + alpha[1L])
    j <- which.max(angles)
    i <- if (j == 1L)
      n
    else j - 1L
    firstandlast <- as.integer(names(angles[c(i, j)]))
    theta0 <- mean(atan2(yy[firstandlast], xx[firstandlast]))
    x0 <- r0 * cos(theta0)
    y0 <- r0 * sin(theta0)
    inc <- diff(pretty(c(0, r0))[1:2])
    srt <- 360 * theta0/(2 * pi)
    coef <- -1
    if (abs(srt) > 90) {
      srt <- srt + 180
      coef <- 1
    }
    len <- 0.025 * r0
    r <- r0
    while (r > 1e-08) {
      x <- r * cos(theta0)
      y <- r * sin(theta0)
      if (len/r < 1) {
        ra <- sqrt(len^2 + r^2)
        thetaa <- theta0 + coef * asin(len/r)
        xa <- ra * cos(thetaa)
        ya <- ra * sin(thetaa)
        segments(xa, ya, x, y)
        text(xa, ya, r0 - r, srt = srt, adj = c(0.5,
                                                1.1), ...)
      }
      r <- r - inc
    }
    segments(x, y, x0, y0)
  }
}
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
########################################################################################################################################

########################################################################################################################################
## Load data
########################################################################################################################################
load('2_analysis_index/1_index_computations/Initial_index_computation_and_parameters_11102023.Rdata')
load('2_analysis_index/2_find_index_groups/Lineages_detected_11102023.Rdata')
########################################################################################################################################



########################################################################################################################################
## 1 - reconstruct ancestral states
########################################################################################################################################
reconstruct_node_states = function(tree, dataset_tips, dataset_with_nodes, min_prop, max_prop, names_seqs){
  # Meta-data with nodes 
  dataset_with_inferred_resonstruction = dataset_with_nodes[,1:4]
  possible_snps = names(dataset_tips[,4:ncol(dataset_tips)])
  
  ## Filter SNP on frequency: only keep SNP that are present in:
  ## >1% of dataset
  ## <99% of dataset
  prevalence = t(apply(dataset_tips[,4:ncol(dataset_tips)], MARGIN = 2, function(x)table(factor(x, levels = c("0", "1")))))
  prevalence_prop = prevalence[,2]/(prevalence[,1]+prevalence[,2])
  prevalence = cbind(prevalence, prevalence_prop)
  possible_snps = names(which(prevalence[,3] >= min_prop & prevalence[,3] <= max_prop))
  
  ## Filter dataset_tips accordingly
  a = which(is.na(match(colnames(dataset_tips), possible_snps)) == F)
  dataset_tips = dataset_tips[,c(1:3, a)]
  
  print(paste0('Going though ', length(possible_snps), ' snps'))
  k=1
  for(i in possible_snps){
    print(paste0(k, ' / ', length(possible_snps), ' snps'))
    
    snp_data = dataset_tips[,which(colnames(dataset_tips) == i)]
    snp_data = as.factor(snp_data)
    names(snp_data) = names_seqs
    
    ## Perform reconstruction for this position
    # rec = phytools::fastAnc(tree = tree, x = snp_data, CI = TRUE, vars = T) ## factAnc lots faster than ace.ML
    rec = ape::ace(phy = tree, x = snp_data, method = "pic")
    rec_all = c(snp_data, rec$ace) ## List of all states: first all tips, then all nodes
    
    ## Find first state, to set it to 0, always
    first_state = rec_all[which(dataset_with_inferred_resonstruction$ID == length(names_seqs) + 1)]
    first_state = round(first_state, digits = 0)
    
    ## Write reconstruction in the big dataset
    if(first_state < 1.5){ ## First state is 0: all good
      dataset_with_inferred_resonstruction = cbind(dataset_with_inferred_resonstruction, rec_all - 1)     
    }
    if(first_state > 1.5){ ## First state is 1: have to change to 0
      dataset_with_inferred_resonstruction = cbind(dataset_with_inferred_resonstruction,  2 - rec_all)
    }
    
    ## Set column name to position name
    colnames(dataset_with_inferred_resonstruction)[which(colnames(dataset_tips) == i) +1] = i
    
    k=k+1
  }
  return(list('dataset_with_inferred_resonstruction' = dataset_with_inferred_resonstruction,
              'snp_prevalence' = prevalence,
              'possible_snps' = possible_snps))
}
########################################################################################################################################

########################################################################################################################################
## Load data
########################################################################################################################################
## Locations
locations = c('World', 'Africa', 'Asia', 'Europe', 'North_America', 'Oceania', 'South_America')

## Timed-tree
tree_sars_all = list(tree_sars_world, 
                     tree_sars_africa, 
                     tree_sars_asia, 
                     tree_sars_europe, 
                     tree_sars_northamerica, 
                     tree_sars_oceania,
                     tree_sars_southamerica)
names(tree_sars_all) = locations

## Names all sequences
names_seqs_all = list(names_seqs_world, 
                      names_seqs_africa, 
                      names_seqs_asia, 
                      names_seqs_europe, 
                      names_seqs_northamerica, 
                      names_seqs_oceania,
                      names_seqs_southamerica)
names(names_seqs_all) = locations

## Collection times of all sequences
times_seqs_all = list(times_seqs_world, 
                      times_seqs_africa, 
                      times_seqs_asia, 
                      times_seqs_europe, 
                      times_seqs_northamerica, 
                      times_seqs_oceania,
                      times_seqs_southamerica)
names(times_seqs_all) = locations

## Correspondences Virus ID // EPI ID (GISAID)
correspondence_ID_EPI_world = read.csv('1_Data_Nextstrain_20230414/World/nextstrain_ncov_gisaid_global_all-time_timetree_virus_ID.csv', header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_asia = read.csv('1_Data_Nextstrain_20230414/Asia/nextstrain_ncov_gisaid_asia_all-time_timetree_virus_ID.csv', header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_africa = read.csv('1_Data_Nextstrain_20230414/Africa/nextstrain_ncov_gisaid_africa_all-time_timetree_virus_ID.csv', header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_europe = read.csv('1_Data_Nextstrain_20230414/Europe/nextstrain_ncov_gisaid_europe_all-time_timetree_virus_ID.csv', header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_northamerica = read.csv('1_Data_Nextstrain_20230414/North_America/nextstrain_ncov_gisaid_north-america_all-time_timetree_virus_ID.csv', header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_oceania = read.csv('1_Data_Nextstrain_20230414/Oceania/nextstrain_ncov_gisaid_oceania_all-time_timetree_virus_ID.csv', header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_southamerica = read.csv('1_Data_Nextstrain_20230414/South_America/nextstrain_ncov_gisaid_south-america_all-time_timetree_virus_ID.csv', header = F, col.names = c("N", "Virus name", "Accession ID"))
correspondence_ID_EPI_all = list(correspondence_ID_EPI_world, 
                                 correspondence_ID_EPI_africa, 
                                 correspondence_ID_EPI_asia, 
                                 correspondence_ID_EPI_europe, 
                                 correspondence_ID_EPI_northamerica, 
                                 correspondence_ID_EPI_oceania,
                                 correspondence_ID_EPI_southamerica)
names(correspondence_ID_EPI_all) = locations

## Datasets with nodes and tips
dataset_with_nodes_all = list(dataset_with_nodes_world, 
                              dataset_with_nodes_africa, 
                              dataset_with_nodes_asia, 
                              dataset_with_nodes_europe, 
                              dataset_with_nodes_northamerica, 
                              dataset_with_nodes_oceania,
                              dataset_with_nodes_southamerica)
names(dataset_with_nodes_all) = locations
########################################################################################################################################

########################################################################################################################################
## Reconstruction
########################################################################################################################################
ORFs = c('E', 'full', 'M', 'N', 'ORF10', 'ORF14', 
         'ORF1a', 'ORF1b', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b',
         'ORF8', 'ORF9b', 'S')

foreach(l = 1:length(locations)) %do% {
  print(paste0('Location: ', locations[l]))
  names_seqs_world_simple = unlist(lapply(names_seqs_all[[l]], function(x){
    tmp = str_split(x, pattern = '\\.')[[1]]
    if(length(tmp) == 1){
      return(tmp)
    }else if(length(tmp) == 2){
      return(tmp[1])
    }else {
      tmp = tmp[-length(tmp)]
      return(paste0(tmp, collapse = '.'))
    }
  }))
  vcf_names = list.files(path = paste0('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/',locations[l], '/', collapse = ''), pattern = '*.vcf')
  
  for(o in 1:length(ORFs)){
    print(paste0('Location: ', locations[l], ' / ORF: ', ORFs[o]))
    data_vcf = read.csv(file = paste0('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/',locations[l], '/', vcf_names[o], collapse = ''), sep = '\t')
    m = match(colnames(data_vcf), correspondence_ID_EPI_all[[l]]$Accession.ID)
    colnames(data_vcf)[which(!is.na(m))] = correspondence_ID_EPI_all[[l]]$Virus.name[m[which(!is.na(m))]]
    
    ## Create dataset with names of each sequence
    dataset_tips = data.frame('ID' = 1:length(names_seqs_all[[l]]),
                              'name_seq' = names_seqs_all[[l]],
                              'time' = times_seqs_all[[l]])
    ## Add AA data to the main dataset
    a = match(names_seqs_world_simple, colnames(data_vcf))
    dataset_tips = cbind(dataset_tips, t(data_vcf[,a]))
    colnames(dataset_tips) = c('ID', 'name_seq', 'time', data_vcf$POS)
    
    ## Reconstruction
    reconstruction = reconstruct_node_states(tree = tree_sars_all[[l]], 
                                             dataset_tips = dataset_tips, 
                                             dataset_with_nodes = dataset_with_nodes_all[[l]], 
                                             min_prop = 0, 
                                             max_prop = 1, 
                                             names_seqs = names_seqs_all[[l]])
    ## Save
    saveRDS(reconstruction, file = paste0('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/', locations[l], '/', locations[l], '_reconstruction_', ORFs[o], '.rds'))
  }
}
########################################################################################################################################







########################################################################################################################################
## 2 - association scores
########################################################################################################################################

########################################################################################################################################
## Load AA data
########################################################################################################################################
## Vcfs
########################################################################################################################################
ORFs = c('E', 'M', 'N', 'ORF10', 'ORF14', 
         'ORF1a', 'ORF1b', 'ORF3a', 'ORF6', 'ORF7a', 'ORF7b',
         'ORF8', 'ORF9b', 'S')

data_vcf_world_E = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_E.vcf", sep = '\t')
data_vcf_world_M = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_M.vcf", sep = '\t')
data_vcf_world_N = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_N.vcf", sep = '\t')
data_vcf_world_ORF10 = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF10.vcf", sep = '\t')
data_vcf_world_ORF14 = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF14.vcf", sep = '\t')
data_vcf_world_ORF1a = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF1a.vcf", sep = '\t')
data_vcf_world_ORF1b = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF1b.vcf", sep = '\t')
data_vcf_world_ORF3a = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF3a.vcf", sep = '\t')
data_vcf_world_ORF6 = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF6.vcf", sep = '\t')
data_vcf_world_ORF7a = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF7a.vcf", sep = '\t')
data_vcf_world_ORF7b = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF7b.vcf", sep = '\t')
data_vcf_world_ORF8 = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF8.vcf", sep = '\t')
data_vcf_world_ORF9b = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_ORF9b.vcf", sep = '\t')
data_vcf_world_S = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_S.vcf", sep = '\t')
data_vcf_world_full = read.csv(file = "2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_full.vcf", sep = '\t')

## Reconstructions
idx_min = which.min(dataset_with_nodes_world$time[which(dataset_with_nodes_world$is.node == 'no')])

reconstruction_world_E = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_E.rds')
reconstruction_world_E$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_E$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_E$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_E$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_E$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_E$possible_snps = colnames(reconstruction_world_E$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_M = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_M.rds')
reconstruction_world_M$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_M$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_M$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_M$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_M$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_M$possible_snps = colnames(reconstruction_world_M$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_N = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_N.rds')
reconstruction_world_N$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_N$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_N$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_N$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_N$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_N$possible_snps = colnames(reconstruction_world_N$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF10 = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF10.rds')
reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF10$possible_snps = colnames(reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF14 = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF14.rds')
reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF14$possible_snps = colnames(reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF1a = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF1a.rds')
reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF1a$possible_snps = colnames(reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF1b = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF1b.rds')
reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF1b$possible_snps = colnames(reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF3a = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF3a.rds')
reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF3a$possible_snps = colnames(reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF6 = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF6.rds')
reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF6$possible_snps = colnames(reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF7a = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF7a.rds')
reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF7a$possible_snps = colnames(reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF7b = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF7b.rds')
reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF7b$possible_snps = colnames(reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF8 = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF8.rds')
reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF8$possible_snps = colnames(reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_ORF9b = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_ORF9b.rds')
reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_ORF9b$possible_snps = colnames(reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon)[-(1:4)]
reconstruction_world_S = readRDS('2_analysis_index/3_snp_association/vcfs_and_AA_reconstructions/World/nextstrain_ncov_gisaid_global_all-time_timetree_reconstruction_S.rds')
reconstruction_world_S$dataset_with_inferred_reconstruction_codon[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_codon)] = reconstruction_world_S$dataset_with_inferred_reconstruction_codon[idx_min,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_codon)]
reconstruction_world_S$dataset_with_inferred_reconstruction_AA[tree_sars_world$Nnode+2,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_AA)] = reconstruction_world_S$dataset_with_inferred_reconstruction_AA[idx_min,5:ncol(reconstruction_world_S$dataset_with_inferred_reconstruction_AA)]
reconstruction_world_S$possible_snps = colnames(reconstruction_world_S$dataset_with_inferred_reconstruction_codon)[-(1:4)]
########################################################################################################################################

########################################################################################################################################
## For each SNP, look at defining mutation of each group
########################################################################################################################################
association_scores_per_group = function(dataset_with_nodes, dataset_with_inferred_reconstruction, tree, 
                                        possible_snps, upstream_window, downstream_window){
  ## Set list to store results
  group_names = levels(as.factor(dataset_with_nodes$groups))
  scores = as.list(rep(NA, length(group_names)-1))
  n_tips = length(tree$tip.label)
  
  ## For each group (except the initial group, which is the root), look at snp association
  for(j in 1:(length(group_names)-1)){ 
    print(j)
    
    ## Find members of the group and MRCA
    members = dataset_with_nodes$ID[which(dataset_with_nodes$groups == group_names[j])]
    mrca = getMRCA(tree, dataset_with_nodes$name_seq[which(dataset_with_nodes$groups == group_names[j] & dataset_with_nodes$is.node == 'no')])
    members = unique(c(members, mrca))
    
    ## Update members list, with strains downstream, within the time window
    time_mrca = dataset_with_nodes$time[which(dataset_with_nodes$ID == mrca)] ## Reference time for window
    # tmp = getDescendants(tree, mrca)
    # tmp = tmp[which(dataset_with_nodes$time[tmp] < time_mrca + downstream_window)]
    # members = unique(c(tmp, members))
    
    ## Update members list, with strains upstream, within the time window
    ancest = Ancestors(x = tree, node = mrca, type = 'all')
    time_ancest = dataset_with_nodes$time[ancest]
    time_ancest[1] = time_mrca ## Always keep the first ancestor node
    tmp = which(time_ancest < time_mrca - upstream_window)
    if(length(tmp) > 0) ancest = ancest[-tmp]  ## Remove nodes that are outside the time window
    if(length(ancest) == 0) ancest = Ancestors(x = tree, node = mrca, type = 'all')[1]
    groups_ancest = dataset_with_nodes$groups[ancest]
    gr = min(as.numeric(as.character(groups_ancest)))
    ancest = ancest[which(groups_ancest == gr)] ## Makes sure we only have the directly ancestral group, not more
    time_ancest = dataset_with_nodes$time[ancest]
    
    ## Branches of interest
    branches = tree$edge
    branches = branches[match(c(ancest, members), branches[,2]),] ## Take all tips, from ancests, members
    tmp = which(is.na(match(branches[,1], c(ancest, members))))
    if(length(tmp) > 0){  branches = branches[-tmp, ]} ## Remove nodes that are not in ancest or members
    
    ## Branches time
    branches_time = branches
    branches_time[,1] = dataset_with_nodes$time[branches[,1]]
    branches_time[,2] = dataset_with_nodes$time[branches[,2]]
    
    ## Branches group (checked: ok)
    branches_group = branches
    branches_group[,1] = dataset_with_nodes$groups[branches[,1]]
    branches_group[,2] = dataset_with_nodes$groups[branches[,2]]
    
    ## Make branch group matrix binary: 1=group of interest, 0=other group (eg ancestral)
    branches_group[which(branches_group == j, arr.ind = T)] = 1
    branches_group[which(branches_group > j, arr.ind = T)] = 0
    
    snps = time_diff = snps_props_within = snps_props_whole = names_possibles_snps = NULL
    
    for(i in 1:length(possible_snps)){
      branches_snp = branches
      branches_snp[,1] = dataset_with_inferred_reconstruction[branches[,1],
                                                              which(colnames(dataset_with_inferred_reconstruction) == possible_snps[i])]
      branches_snp[,2] = dataset_with_inferred_reconstruction[branches[,2],
                                                              which(colnames(dataset_with_inferred_reconstruction) == possible_snps[i])]
      
      ancestral_state = branches_snp[which.min(branches[,1]),1]
      
      k=1
      ancestral_state = branches_snp[which(branches[,1] == rev(ancest)[k]),1]
      k=2
      while(str_detect(ancestral_state, pattern = 'n|X') == T & k <= length(ancest)-1){
        ancestral_state = branches_snp[which(branches[,1] == rev(ancest)[k]),1]
        k=k+1
      }
      if(str_detect(ancestral_state, pattern = 'n|X') == T){
        ancestral_state = branches_snp[which(branches[,1] == mrca),1][1]
      }
      
      branches_snp[str_detect(branches_snp, pattern = 'n|X')] = ancestral_state
      
      branches_snp_bin = (branches_snp!=ancestral_state)
      
      tmp = which(branches_group[,1] == 1 & branches_group[,2] == 1)
      Px = as.numeric(branches_group[tmp,])
      Sx = as.numeric(branches_snp_bin[tmp,])
      score3 <- (sum((1 - Px)*(1 - Sx), na.rm=TRUE) + sum(Px*Sx, na.rm=TRUE)) / length(Px)
      
      scores[[j]][i] = score3
      
      t = table(branches_snp)
      t = t[which(names(t) != ancestral_state)]
      t = t[which.max(t)]
      names_possibles_snps = c(names_possibles_snps, paste0(ancestral_state, '|', possible_snps[i], '|', names(t)))
    }
    scores[[j]] = scores[[j]]#- median(scores[[j]])
    names(scores[[j]]) = names_possibles_snps
  }
  return(scores)
}
## Time window
upstream_window = 2 # years (not going to the root)
downstream_window = 10 # years (considering all the group)

scores_world_E_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_E$dataset_with_inferred_reconstruction_codon, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_E$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_M_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_M$dataset_with_inferred_reconstruction_codon, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_M$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_N_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_N$dataset_with_inferred_reconstruction_codon, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_N$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF10_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF10$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF14_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF14$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF1a_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF1a$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF1b_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF1b$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF3a_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF3a$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF6_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                        dataset_with_inferred_reconstruction = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_codon, 
                                                        tree = tree_sars_world, 
                                                        possible_snps = reconstruction_world_ORF6$possible_snps, 
                                                        upstream_window = upstream_window, 
                                                        downstream_window = downstream_window)
scores_world_ORF7a_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF7a$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF7b_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF7b$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_ORF8_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                        dataset_with_inferred_reconstruction = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_codon, 
                                                        tree = tree_sars_world, 
                                                        possible_snps = reconstruction_world_ORF8$possible_snps, 
                                                        upstream_window = upstream_window, 
                                                        downstream_window = downstream_window)
scores_world_ORF9b_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                         dataset_with_inferred_reconstruction = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_codon, 
                                                         tree = tree_sars_world, 
                                                         possible_snps = reconstruction_world_ORF9b$possible_snps, 
                                                         upstream_window = upstream_window, 
                                                         downstream_window = downstream_window)
scores_world_S_codons = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_S$dataset_with_inferred_reconstruction_codon, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_S$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)

scores_world_E_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                 dataset_with_inferred_reconstruction = reconstruction_world_E$dataset_with_inferred_reconstruction_AA, 
                                                 tree = tree_sars_world, 
                                                 possible_snps = reconstruction_world_E$possible_snps, 
                                                 upstream_window = upstream_window, 
                                                 downstream_window = downstream_window)
scores_world_M_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                 dataset_with_inferred_reconstruction = reconstruction_world_M$dataset_with_inferred_reconstruction_AA, 
                                                 tree = tree_sars_world, 
                                                 possible_snps = reconstruction_world_M$possible_snps, 
                                                 upstream_window = upstream_window, 
                                                 downstream_window = downstream_window)
scores_world_N_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                 dataset_with_inferred_reconstruction = reconstruction_world_N$dataset_with_inferred_reconstruction_AA, 
                                                 tree = tree_sars_world, 
                                                 possible_snps = reconstruction_world_N$possible_snps, 
                                                 upstream_window = upstream_window, 
                                                 downstream_window = downstream_window)
scores_world_ORF10_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF10$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF10$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF14_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF14$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF14$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF1a_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF1a$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF1a$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF1b_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF1b$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF1b$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF3a_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF3a$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF3a$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF6_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                    dataset_with_inferred_reconstruction = reconstruction_world_ORF6$dataset_with_inferred_reconstruction_AA, 
                                                    tree = tree_sars_world, 
                                                    possible_snps = reconstruction_world_ORF6$possible_snps, 
                                                    upstream_window = upstream_window, 
                                                    downstream_window = downstream_window)
scores_world_ORF7a_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF7a$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF7a$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF7b_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF7b$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF7b$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_ORF8_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                    dataset_with_inferred_reconstruction = reconstruction_world_ORF8$dataset_with_inferred_reconstruction_AA, 
                                                    tree = tree_sars_world, 
                                                    possible_snps = reconstruction_world_ORF8$possible_snps, 
                                                    upstream_window = upstream_window, 
                                                    downstream_window = downstream_window)
scores_world_ORF9b_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                     dataset_with_inferred_reconstruction = reconstruction_world_ORF9b$dataset_with_inferred_reconstruction_AA, 
                                                     tree = tree_sars_world, 
                                                     possible_snps = reconstruction_world_ORF9b$possible_snps, 
                                                     upstream_window = upstream_window, 
                                                     downstream_window = downstream_window)
scores_world_S_AA = association_scores_per_group(dataset_with_nodes = dataset_with_nodes_world,
                                                 dataset_with_inferred_reconstruction = reconstruction_world_S$dataset_with_inferred_reconstruction_AA, 
                                                 tree = tree_sars_world, 
                                                 possible_snps = reconstruction_world_S$possible_snps, 
                                                 upstream_window = upstream_window, 
                                                 downstream_window = downstream_window)
save.image('Raw_scores_detection_17122023.Rdata')
########################################################################################################################################

########################################################################################################################################
## Find significant snps
########################################################################################################################################
edge_lineage_tree = split_world$lineage_tree$edge
edge_lineage_tree[,1] = split_world$tip_and_nodes_groups[match(edge_lineage_tree[,1],names(split_world$tip_and_nodes_groups))]
edge_lineage_tree[,2] = split_world$tip_and_nodes_groups[match(edge_lineage_tree[,2],names(split_world$tip_and_nodes_groups))]
edge_lineage_tree_snps = as.list(edge_lineage_tree)

## Combine all codons
scores_world_sig_all_codons = scores_world_S_codons
sig_threshold = 0.8
for(i in 1:length(scores_world_S_codons)){
  tmp = as.numeric(edge_lineage_tree[which(edge_lineage_tree[,2] == i),])
  if(tmp[1] == length(scores_world_S_codons)+1) {
    ans = NULL
  }else{
    ans = c(paste0('E:', names(which(scores_world_E_codons[[tmp[1]]] > sig_threshold))),
            paste0('M:', names(which(scores_world_M_codons[[tmp[1]]] > sig_threshold))),
            paste0('N:', names(which(scores_world_N_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF10:', names(which(scores_world_ORF10_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF14:', names(which(scores_world_ORF14_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF1a:', names(which(scores_world_ORF1a_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF1b:', names(which(scores_world_ORF1b_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF3a:', names(which(scores_world_ORF3a_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF6:', names(which(scores_world_ORF6_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF7a:', names(which(scores_world_ORF7a_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF7b:', names(which(scores_world_ORF7b_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF8:', names(which(scores_world_ORF8_codons[[tmp[1]]] > sig_threshold))),
            paste0('ORF9b:', names(which(scores_world_ORF9b_codons[[tmp[1]]] > sig_threshold))),
            paste0('S:', names(which(scores_world_S_codons[[tmp[1]]] > sig_threshold))))
  }
  des = c(paste0('E:', names(which(scores_world_E_codons[[tmp[2]]] > sig_threshold))),
          paste0('M:', names(which(scores_world_M_codons[[tmp[2]]] > sig_threshold))),
          paste0('N:', names(which(scores_world_N_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF10:', names(which(scores_world_ORF10_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF14:', names(which(scores_world_ORF14_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF1a:', names(which(scores_world_ORF1a_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF1b:', names(which(scores_world_ORF1b_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF3a:', names(which(scores_world_ORF3a_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF6:', names(which(scores_world_ORF6_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF7a:', names(which(scores_world_ORF7a_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF7b:', names(which(scores_world_ORF7b_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF8:', names(which(scores_world_ORF8_codons[[tmp[2]]] > sig_threshold))),
          paste0('ORF9b:', names(which(scores_world_ORF9b_codons[[tmp[2]]] > sig_threshold))),
          paste0('S:', names(which(scores_world_S_codons[[tmp[2]]] > sig_threshold))))
  res = des[which(is.na(match(des, ans)))]
  res = res[which(is.na(match(res, c("E:", "M:", "N:", 'ORF10:', 'ORF14:',
                                     'ORF1a:', 'ORF1b:', 'ORF3a:', 'ORF6:', 
                                     'ORF7a:', 'ORF7b:', 'ORF8:', 
                                     'ORF9b:','S:'))))]
  scores_world_sig_all_codons[[i]] = res
}

## Combine all AA
scores_world_sig_all_AA = scores_world_S_AA
sig_threshold = 0.8
for(i in 1:length(scores_world_S_AA)){
  tmp = as.numeric(edge_lineage_tree[which(edge_lineage_tree[,2] == i),])
  if(tmp[1] == length(scores_world_S_AA)+1) {
    ans = NULL
  }else{
    ans = c(paste0('E:', names(which(scores_world_E_AA[[tmp[1]]] > sig_threshold))),
            paste0('M:', names(which(scores_world_M_AA[[tmp[1]]] > sig_threshold))),
            paste0('N:', names(which(scores_world_N_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF10:', names(which(scores_world_ORF10_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF14:', names(which(scores_world_ORF14_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF1a:', names(which(scores_world_ORF1a_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF1b:', names(which(scores_world_ORF1b_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF3a:', names(which(scores_world_ORF3a_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF6:', names(which(scores_world_ORF6_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF7a:', names(which(scores_world_ORF7a_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF7b:', names(which(scores_world_ORF7b_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF8:', names(which(scores_world_ORF8_AA[[tmp[1]]] > sig_threshold))),
            paste0('ORF9b:', names(which(scores_world_ORF9b_AA[[tmp[1]]] > sig_threshold))),
            paste0('S:', names(which(scores_world_S_AA[[tmp[1]]] > sig_threshold))))
  }
  des = c(paste0('E:', names(which(scores_world_E_AA[[tmp[2]]] > sig_threshold))),
          paste0('M:', names(which(scores_world_M_AA[[tmp[2]]] > sig_threshold))),
          paste0('N:', names(which(scores_world_N_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF10:', names(which(scores_world_ORF10_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF14:', names(which(scores_world_ORF14_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF1a:', names(which(scores_world_ORF1a_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF1b:', names(which(scores_world_ORF1b_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF3a:', names(which(scores_world_ORF3a_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF6:', names(which(scores_world_ORF6_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF7a:', names(which(scores_world_ORF7a_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF7b:', names(which(scores_world_ORF7b_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF8:', names(which(scores_world_ORF8_AA[[tmp[2]]] > sig_threshold))),
          paste0('ORF9b:', names(which(scores_world_ORF9b_AA[[tmp[2]]] > sig_threshold))),
          paste0('S:', names(which(scores_world_S_AA[[tmp[2]]] > sig_threshold))))
  res = des[which(is.na(match(des, ans)))]
  res = res[which(is.na(match(res, c("E:", "M:", "N:", 'ORF10:', 'ORF14:',
                                     'ORF1a:', 'ORF1b:', 'ORF3a:', 'ORF6:', 
                                     'ORF7a:', 'ORF7b:', 'ORF8:', 
                                     'ORF9b:','S:'))))]
  scores_world_sig_all_AA[[i]] = res
}
########################################################################################################################################

########################################################################################################################################
## Full genome, substitutions (NS, ie AA change) with signal
########################################################################################################################################
genes = read.csv('1_refseq/Position_ORFs_nextstrain.csv')
genes = rbind(genes, 
              c(15, NA, NA, 'region', 21563+(333-1)*3, 21563+(527-1)*3+2, NA, '+', NA, 'region_name "RBD"', "RBD"))
genes_to_plot = genes
genes_to_plot = genes_to_plot[-which(genes_to_plot$gene == 'ORF3a' | genes_to_plot$gene == 'ORF10' | genes_to_plot$gene == 'ORF14'| genes_to_plot$gene == 'ORF6' | genes_to_plot$gene == 'ORF7a'|
                                       genes_to_plot$gene == 'ORF7b' | genes_to_plot$gene == 'ORF8'| genes_to_plot$gene == 'ORF9b'),]
colors = MetBrewer::met.brewer(name="Derain", n=15)

scores_world_full_AA = NULL
for(i in 1:(nrow(genes)-1)){
  if(genes$gene[i] == 'E') {
    tmp = apply(do.call(rbind, scores_world_E_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'M') {
    tmp = apply(do.call(rbind, scores_world_M_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'N') {
    tmp = apply(do.call(rbind, scores_world_N_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF10') {
    tmp = apply(do.call(rbind, scores_world_ORF10_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF14') {
    tmp = apply(do.call(rbind, scores_world_ORF14_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF1a') {
    tmp = apply(do.call(rbind, scores_world_ORF1a_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF1b') {
    tmp = apply(do.call(rbind, scores_world_ORF1b_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF3a') {
    tmp = apply(do.call(rbind, scores_world_ORF3a_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF6') {
    tmp = apply(do.call(rbind, scores_world_ORF6_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF7a') {
    tmp = apply(do.call(rbind, scores_world_ORF7a_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF7b') {
    tmp = apply(do.call(rbind, scores_world_ORF7b_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF8') {
    tmp = apply(do.call(rbind, scores_world_ORF8_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'ORF9b') {
    tmp = apply(do.call(rbind, scores_world_ORF9b_AA), MAR = 2, max)
  } else if (genes$gene[i] == 'S') {
    tmp = apply(do.call(rbind, scores_world_S_AA), MAR = 2, max)
  } 
  names_tmp = as.numeric(unlist(lapply(names(tmp), function(x)str_split(x, '\\|')[[1]][2])))
  names_tmp = names_tmp*3-1.5 + as.numeric(genes$start[i]) - 1
  names(tmp) = names_tmp
  scores_world_full_AA = c(scores_world_full_AA, tmp)
}

plot_scores_genome_density_AA = function(){
  par(oma = c(0,0,0,0), mar = c(2,2,0,0), mgp = c(0,0.1,0), family = 'Arial',
      cex.axis=0.3, cex.lab=0.3, cex.main=0.4, cex.sub=0.3)
  dens = density(as.numeric(names(scores_world_full_AA)), weights = scores_world_full_AA/sum(scores_world_full_AA), 
                 bw = 50, n = 3000)
  plot(NULL, xlim = c(0, 30000), ylim = c(0, max(dens$y)), bty = 'n', yaxt = 'n',
       ylab = '')
  axis(2, las = 2)
  for(i in 1:nrow(genes)){
    if(!is.na(match(genes$gene[i], genes_to_plot$gene))){
      polygon(x = c(c(genes$start[i], genes$end[i]), rev(c(genes$start[i], genes$end[i]))), y = c(0, 0, 1, 1), 
              border = F, col = adjustcolor(colors[i], alpha.f = 0.25))
    }
  }
  polygon(x = c(dens$x, rev(dens$x)), 
          y = c(dens$y, rep(0, length(dens$y))), 
          border = F, col = 'grey30')
}
ggdraw(plot_scores_genome_density_AA)
ggsave(filename = 'Plot_score_AA_genome_density_20231109.pdf', 
       plot = ggdraw(plot_scores_genome_density_AA), 
       device = 'pdf', scale = 1,
       width = 10, height = 5, units = 'cm')
############################################

########################################################################################################################################
## Save list lineage-defining mutations
########################################################################################################################################
df = data.frame('Lineage' = 1:13, 
                'Mutations' = NA)
df$Mutations = unlist(lapply(scores_world_sig_all_AA, function(x)paste0(x, collapse = '; ')))

write.csv(df, 'List_AA_defining_lineages.csv')
########################################################################################################################################

