# clear R's brain
rm(list = ls())


# packages

library(TDbook)
library(ggtree)
library(phangorn)
library(tidyverse)
library(phytools)



# info here: https://yulab-smu.top/treedata-book/chapter7.html

setwd("/Users/habich/Documents/Projects/AH123_accessory_effector_evolution/b_mirror_trees/")

eff <- 'PA14-43100_PA14'

tree_wg <- read.newick(paste0(eff,'/', eff, '_core_gene_alignment_filtered.aln.treefile'))
tree_ig <- read.newick(paste0(eff,'/',eff,'_msa.fasta.treefile'))

tree_ig_new_tip_label <- gsub(paste0(eff,'_'), '', tree_ig$tip.label)
tree_ig$tip.label <- tree_ig_new_tip_label

tips_wg <- tree_wg$tip.label
tips_ig <- tree_ig$tip.label

#creation of the association matrix:
association <- cbind(sort(tree_wg$tip.label), sort(tree_ig$tip.label))

tree_wg2 <- midpoint(tree_wg)
tree_ig2 <- midpoint(tree_ig)

tree_co_phytools <- cophylo(tree_wg2, tree_ig2,
                            assoc = association,
                            size = 0.05)

plot(tree_co_phytools, scale.bar=c(0.01,0.01), fsize=0.2,
     link.lty='solid', link.lwd=0.05)
