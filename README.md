Code for the paper ‘Learning the fitness dynamics of pathogens from
phylogenies’
================
Noémie Lefrancq
January 2 2024

- <a href="#content-of-this-repo" id="toc-content-of-this-repo">Content of
  this repo</a>
- <a href="#example-on-sars-cov-2" id="toc-example-on-sars-cov-2">Example
  on SARS-CoV-2</a>
  - <a href="#load-codes-and-data" id="toc-load-codes-and-data">Load codes
    and data</a>
  - <a href="#compute-index" id="toc-compute-index">Compute index</a>
  - <a href="#plot-tree--index-below-with-colors-from-nextstrain-clades"
    id="toc-plot-tree--index-below-with-colors-from-nextstrain-clades">Plot
    tree &amp; index below, with colors from NextStrain clades</a>
  - <a href="#find-clades-based-on-index-dynamics"
    id="toc-find-clades-based-on-index-dynamics">Find clades based on index
    dynamics</a>
  - <a href="#plot-tree--index-below-with-colors-from-index-defined-groups"
    id="toc-plot-tree--index-below-with-colors-from-index-defined-groups">Plot
    tree &amp; index below, with colors from index-defined groups</a>
  - <a href="#compare-nextstrain-groups-and-groups-called-with-the-index"
    id="toc-compare-nextstrain-groups-and-groups-called-with-the-index">Compare
    NextStrain groups and groups called with the index</a>

# Content of this repo

In the repo you will find:

1.  In `1_Data`, the data used in the paper.
2.  In `2_Functions`, the codes that are behind the analysis:
    - fitness index dynamics
    - lineages detection
    - lineage fitness estimation
    - lineage-defining mutations
3.  In `3_Analysis_per_pathogen`, the analysis per pathogen

*Below is a brief example of what can be archived on SARS-CoV-2, using
the codes in the folder `2_Codes`.*

# Example on SARS-CoV-2

## Load codes and data

#### Load index functions

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

``` r
## source XXX
```

#### Load toy data

Loading SARS-CoV-2 tree, in which all the tip name include: collection
time, location and Pango lineage

``` r
library(ape) ; library(stringr)
tree_sars_cov2 = read.nexus('~/Documents/THD/Index_codes/Data/nextstrain_20220913_ncov_gisaid_global_all-time_timetree_tipdates_tipclades.nexus')
## Make sure the tree is binary, and ladderized
tree_sars_cov2 = collapse.singles(ladderize(multi2di(tree_sars_cov2)))
## Names all sequences
names_seqs = tree_sars_cov2$tip.label
## Collection times of all sequences
times_seqs = as.numeric(sapply(names_seqs, function(x)tail(str_split(x, pattern = '/')[[1]],2)[1]))
## Nextstrain clades of all sequences
clades_seqs = sapply(names_seqs, function(x)tail(str_split(x, pattern = '/')[[1]],1))
```

#### Index parameters

``` r
## Length genome 
genome_length = 29903 # reference nextstrain https://www.ncbi.nlm.nih.gov/nuccore/MN908947
## Mutation rate 
mutation_rate = 8.1e-4 # mutation rate used by nextstrain https://github.com/nextstrain/ncov
## Parameters for the index
timescale = 0.15 ## Timescale
## Window of time on which to search for samples in the population
wind = 15 #days
wind = wind/365
```

## Compute index

#### Compute pairwise distance matrix

Compute distance between each pair of sequences and nodes in the tree

``` r
genetic_distance_mat = dist.nodes.with.names(tree_sars_cov2)
```

Get the time of each node

``` r
nroot = length(tree_sars_cov2$tip.label) + 1 ## Root number
distance_to_root = genetic_distance_mat[nroot,]
root_height = times_seqs[which(names_seqs == names(distance_to_root[1]))] - distance_to_root[1]
nodes_height = root_height + distance_to_root[length(names_seqs)+(1:(length(names_seqs)-1))]
```

#### Preparation data tips and nodes

``` r
# Meta-data with nodes 
dataset_with_nodes = data.frame('ID' = c(1:length(names_seqs), length(names_seqs)+(1:(length(names_seqs)-1))),
                                'name_seq' = c(names_seqs, length(names_seqs)+(1:(length(names_seqs)-1))),
                                'time' = c(times_seqs, nodes_height),
                                'is.node' = c(rep('no', length(names_seqs)), rep('yes', (length(names_seqs)-1))),
                                'clade' = c(clades_seqs, rep(NA, length(names_seqs)-1)))
```

#### Compute index of every tip and node

``` r
dataset_with_nodes$index = Index::compute.index(time_distance_mat = genetic_distance_mat, 
                                                timed_tree = tree_sars_cov2, 
                                                time_window = wind,
                                                metadata = dataset_with_nodes, 
                                                mutation_rate = mutation_rate,
                                                timescale = timescale,
                                                genome_length = genome_length)
```

## Plot tree & index below, with colors from NextStrain clades

First, generate the color key, based on the Nextstrain clade of each
sequence.

``` r
library(MetBrewer)
colors_clade = met.brewer(name="Cross", n=length(levels(as.factor(dataset_with_nodes$clade))), type="continuous")

dataset_with_nodes$color = as.factor(dataset_with_nodes$clade)
clade_labels = levels(dataset_with_nodes$color)
levels(dataset_with_nodes$color) = colors_clade
dataset_with_nodes$color = as.character(dataset_with_nodes$color)
```

Then plot the tree and index:

``` r
par(mfrow = c(2,1), oma = c(0,0,0,0), mar = c(4,4,0,0))

## Tree
plot(tree_sars_cov2, show.tip.label = FALSE, edge.color = 'grey', edge.width = 0.25)
tiplabels(pch = 16, col = dataset_with_nodes$color, cex = 0.25)
nodelabels(pch = 16, col = 'grey', cex = 0.25)
axisPhylo(side = 1, root.time = root_height, backward = F)

## Index
plot(dataset_with_nodes$time, 
     dataset_with_nodes$index, 
     col = adjustcolor(dataset_with_nodes$color, alpha.f = 1),
     bty = 'n', xlim = c(2019.95, 2022.7), cex = 0.5,
     pch = 16, bty = 'n', ylim = c(0, 1), 
     main = paste0(''), 
     ylab = 'Diversity index', xlab = 'Time (years)', yaxt = 'n')
axis(2, las = 2)

# Color key
legend('topright', 
       legend = clade_labels,
       fill = colors_clade, border = colors_clade,
       cex = 0.5, bty = 'n', ncol = 5)
```

## Find clades based on index dynamics

#### Compute index of every tip and node

Find splits (can be quite long)

``` r
potential_splits = find.groups.by.index.dynamics(timed_tree = tree_sars_cov2, 
                                                 time_distance_mat = genetic_distance_mat, 
                                                 metadata = dataset_with_nodes,
                                                 mutation_rate = mutation_rate, 
                                                 genome_length = genome_length, 
                                                 min_offspring_per_start_nodes = 100,
                                                 time_window = 1, 
                                                 min_number_clade = 5, 
                                                 p_value_threshold = 0.01, 
                                                 error_threshold = 0.05,
                                                 show_progress = T)
```

Optimize the number of groups: set the minimum number of sequences per
group to 20, with a minimum frequency of 10%.

``` r
split = merge.groups(timed_tree = tree_sars_cov2, metadata = dataset_with_nodes, 
                     initial_splits = potential_splits, 
                     group_count_threshold = 20, group_freq_threshold = 0.05)
```

Label sequences with these new groups, and assign a color to each of
them.

``` r
library(MetBrewer)
## Label sequences with new groups
dataset_with_nodes$groups = as.factor(split$groups)
name_groups = levels(dataset_with_nodes$groups)
n_groups <- length(name_groups)

## Choose color palette
colors_groups = met.brewer(name="Cross", n=n_groups, type="continuous")

## Color each group
dataset_with_nodes$group_color = dataset_with_nodes$groups
levels(dataset_with_nodes$group_color) = colors_groups
dataset_with_nodes$group_color = as.character(dataset_with_nodes$group_color)
```

## Plot tree & index below, with colors from index-defined groups

Plot the tree and index colored with the new groups:

``` r
par(mfrow = c(2,1), oma = c(0,0,0,0), mar = c(4,4,0,0))

## Tree
plot(tree_sars_cov2, show.tip.label = FALSE, 
     edge.color = 'grey', edge.width = 0.25)
tiplabels(pch = 16, col = dataset_with_nodes$group_color, cex = 0.5)
axisPhylo(side = 1, root.time = root_height, backward = F)

## Index colored by group
plot(dataset_with_nodes$time, 
     dataset_with_nodes$index, 
     col = adjustcolor(dataset_with_nodes$group_color, alpha.f = 1),
     bty = 'n', xlim = c(2019.95, 2022.7), cex = 0.5,
     pch = 16, bty = 'n', ylim = c(0, 1), 
     main = paste0(''), 
     ylab = 'Diversity index', xlab = 'Time (years)', yaxt = 'n')
axis(2, las = 2)
# Color key
legend('topright', 
       legend = name_groups,
       fill = colors_groups, border = colors_groups,
       cex = 0.5, bty = 'n', ncol = 5)
```

Plot the relationships between groups:

``` r
suppressMessages(library(ggtree, quiet = T))
library(ggplot2, quiet = T)
order_colors = order(as.numeric(split$tip_and_nodes_groups))
p = ggtree(split$lineage_tree, aes(col = split$tip_and_nodes_groups), 
           layout="roundrect", size=1.5) +
  geom_point(size=5, alpha=1, aes(col = split$tip_and_nodes_groups)) +
  scale_color_manual(values = colors_groups[match(split$tip_and_nodes_groups[order_colors], name_groups)],
                     breaks = split$tip_and_nodes_groups[order_colors], 
                     na.value = 'white', name = 'Groups',
                     guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=3))
p    
```

## Compare NextStrain groups and groups called with the index

Plot SARS-CoV-2 trees colored with each set of groups next to each
other:

``` r
par(mfrow = c(1,2), oma = c(0,0,0,0), mar = c(4,1,1,0))

## Tree with NextStrain clades
plot(tree_sars_cov2, show.tip.label = FALSE, edge.color = 'grey', edge.width = 0.25,
     main = 'NextStrain clades')
tiplabels(pch = 16, col = dataset_with_nodes$color, cex = 0.5)
axisPhylo(side = 1, root.time = root_height, backward = F)

## Tree with index-defined groups
plot(tree_sars_cov2, show.tip.label = FALSE, edge.color = 'grey', edge.width = 0.25,
     main = 'Index automatic groups',
     direction = "leftwards")
tiplabels(pch = 16, col = dataset_with_nodes$group_color, cex = 0.5)
axisPhylo(side = 1, root.time = root_height, backward = F)
```
