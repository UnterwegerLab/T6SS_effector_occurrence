# clear R's brain
rm(list = ls())


# packages

library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)



# info here: https://www.youtube.com/watch?v=ywHVb0Q-qsM

setwd("/Users/habich/Documents/Projects/AH69_ES_large_scale/downstream_analysis_paper/g_rarefaction_curve/")

g2_md_pre <- read.csv('../forEff_effector_sets.csv', header = T)

g2_md_pre$number_of_eff <- NULL

# define dataset/phylogroup possibilities: A, B, all
g2_ds <- 'all'

# define maximum number of genomes in phylogroup for maximum number of genomes per curve, when comparing phylogroup A with phylgroup B

#md_temp <- md[md_pre$Phylogroup == 'C1',]
#max_B <- length(c(md_temp$Sequence_ID))
g2_md <- unite(g2_md_pre, ES, 2:27, sep = '-')

# subset dataset
# kommentiere die nächste Zeile aus, wenn ds = all
#md <- md[md$phylogroup == ds,]
g2_md$phylogroup <- NULL

#prepare dataset: unite all effectors into one string = effector set


# number of collection curves (the higher the better)
g2_cc <- 10000

#number of genomes per curve (can also be max_B)
g2_ng <- 1912


#___________________________________

# do collection curve several times

g2_output_df <- data.frame(c(1:g2_ng))
colnames(g2_output_df) <- c('number_of_genomes')

for (g2_i in 1:g2_cc) {
  #g2_i <- 1
  # choose ng random genomes
  g2_rg <- sample(g2_md$Sequence_ID, g2_ng, replace = F)
  g2_rg_df <- g2_md[(g2_md$Sequence_ID %in% g2_rg),]
  g2_rg_df$Sequence_ID <- NULL
  
  # generate random order and then arrange by it
  g2_order_v <- sample(1:g2_ng, g2_ng, replace = F)
  g2_rg_df$draw <- g2_order_v
  g2_rg_df <- g2_rg_df %>% arrange(draw)
  g2_collect_df <- data.frame()
  
  #colnames(collect_df) <- c('number_of_genomes', paste('number_of_different_effectors_iteration_', i, sep=''))
  g2_sub_v <- c(1)
  for (g2_j in 2:(g2_ng+1)) {
    #j <- 2
    g2_sub_df <- g2_rg_df[(g2_rg_df$draw %in% g2_sub_v),]
    
    #count how many different effector sets
    g2_sub_df <- melt(g2_sub_df, id = 'draw')
    g2_sub_df <- g2_sub_df[(g2_sub_df$value != 'X'),]
    g2_diff_eff <- length(unique(g2_sub_df$value))

    # add to output dataframe
    g2_new_df <- data.frame(g2_j - 1, g2_diff_eff)
    colnames(g2_new_df) <- c('number_of_genomes', paste('number_of_different_effector_sets_iteration_', g2_i, sep=''))
    g2_collect_df <- rbind(g2_collect_df, g2_new_df)
    
    g2_sub_v <- c(g2_sub_v, g2_j)
  }

  g2_output_df <- dplyr::full_join(g2_output_df, g2_collect_df, by = 'number_of_genomes')

}

g2_output_df$rarefaction <- rowSums(g2_output_df[,2:(g2_cc+1)]) / g2_cc

# plot

g2_melt <- g2_output_df
g2_melt$rarefaction <- NULL
g2_melt <- melt(g2_melt, id = 'number_of_genomes')
colnames(g2_melt) <- c('number_of_genomes', 'number_of_iteration', 'number_of_effector_sets')

ggplot(g2_melt, aes(x = number_of_genomes, y = number_of_effector_sets)) +
  geom_point(shape = 20, color = 'gray90', alpha = 0.1) +
  stat_summary(
    geom = "line",
    fun = "mean",
    col = "black") +
  theme_bw() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 400), breaks = c(0, 100, 200, 300, 400)) +
  scale_x_continuous(expand = c(0,0), limits = c(1, 1915), breaks = c(0, 500, 1000, 1500, 1912)) +
  ylab('Number of different effector sets') +
  xlab('Size of genome subsets') +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_blank()) +
  theme(axis.line = element_line(color = "black", linewidth = 0.25, linetype = "solid"))
ggsave('g2_paper_figure.png', width = 8, height = 8, unit = c("cm"))
ggsave('g2_paper_figure.tiff', device = 'tiff', dpi = 1000)




