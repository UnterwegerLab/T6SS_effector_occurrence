# clear R's brain
rm(list = ls())


library(ggplot2)
library(UpSetR)
library(dplyr)
library(reshape2)
library(packcircles)


setwd("/Users/habich/Documents/Projects/AH69_ES_large_scale/downstream_analysis_paper/f_effector_sets/")

md <- read.csv('../forEff_effector_sets.csv', header = T)

loci <- c('PA1844', 'PA2684', 'PA2702', 'PA2774', 'PA3484', 'PA4163',
          'PA0093', 'PA0099', 
          'PA1510', 'PA1863', 'PA2066', 'PA3290', 'PA3907', 'PA4922', 'PA5089',
          'PA0260', 'PALA36_00282', 'PA0262', 'PA0822', 'PA3487', 'PA14_43100', 'PAKAF_03765', 'PALA52_00191', 'PA5265', 
          'PA0256', 'PA2374',
          'PA14_33970')

effectors <- c('X',
               'PA0256', 'PA2066', 'tplE', 'tse1', 'modA', 'tseF', 'tse5', 'tse2', 'tse4', 'tle1', 'tse3', 'tseT', 'tse8', 'azu', 'pldB',
               'tle3', 'tle4b', 'vgrG2b', 'tseV', 'pldA', 'rhsP2', 'tle2', 'tspE1a',
               'tse6', 'tas1', 'tse6a',
               'tse7', 'tse7a', 'tse7b', 'tse7c', 'tse7d', 'tse7e',
               'tspE1b', 'tspE1c',
               'tepB', 'tepBa', 'tepBb', 'tepBc')

mycols <- c('white',
            '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', '#808080', 
            '#c5c5c5', '#c5c5c5', '#c5c5c5', '#c5c5c5', '#c5c5c5', '#c5c5c5', '#c5c5c5', '#c5c5c5',
            "#c6ddc2", "#82ba6a", "#4c7d4e",
            "#dce7e9", "#a9d6ed", "#89aecf", "#617c94", "#43585a", "#412f6e",
            "#da7a60", "#844633",
            "#f6f1bd", "#f5eb5e", "#c4b353", "#6f6432")

#split up dataset into H1 and H2

H1 <- c('PA1844', 'PA2684', 'PA2702', 'PA2774', 'PA3484', 'PA4163', 'PA0093', 'PA0099')
H2 <- c('PA1510', 'PA1863', 'PA2066', 'PA3290', 'PA3907', 'PA4922', 'PA5089', 'PA0260', 'PALA36_00282', 'PA0262', 'PA0822', 'PA3487', 'PA14_43100', 'PAKAF_03765', 'PALA52_00191', 'PA5265')
H3 <- c('PA0256', 'PA2374', 'PA14_33970')


# ------------- H1 ---------------
md_H1 <- subset(md, select = H1)
#### distinct effector sets
distinct_ES <- md_H1

# identify distinct ES and count how often they occur
distinct_ES <- distinct_ES %>% group_by(across()) %>% count

# give name according to prevalence
distinct_ES <- distinct_ES[order(desc(distinct_ES$n)),]
distinct_ES$ES <- paste('H1_ES', as.character(1:nrow(distinct_ES)), sep = '')
distinct_ES_H1 <- distinct_ES
write.csv(distinct_ES, 'H1_distinct_ES.csv', row.names = F)


# figure: bubble chart
# info: https://cran.r-project.org/web/packages/packcircles/vignettes/progressive_packing.html
md_bubble <- data.frame(distinct_ES$ES, distinct_ES$n)
colnames(md_bubble) <- c('ES', 'n')

bubble_v <- c(md_bubble$n)
packing <- circleProgressiveLayout(bubble_v)
md_bubble_fig <- circleLayoutVertices(packing, npoints = length(bubble_v))

ggplot(md_bubble_fig) +
  geom_polygon(aes(x,y,group = id)) +
  coord_equal() +
  theme_bw() + 
  scale_y_reverse() +
  theme(panel.grid = element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        panel.border = element_blank())
ggsave('H1_bubbles.pdf',  width = 10, height = 10, units = c("cm"))

# to get different layouts
ncircles <- length(bubble_v) 
packings <- lapply(1,
  function(i) { 
    x <- sample(bubble_v, ncircles) 
    circleProgressiveLayout(x) 
  })

packings <- do.call(rbind, packings)
npts <- length(bubble_v)
dat_gg <- circleLayoutVertices(packings, npoints = npts) 

ggplot(data = dat_gg, aes(x, y)) + 
  geom_polygon(aes(group = id),
               fill = 'maroon', color = 'white') +
  coord_equal() +
  theme_bw() +
  scale_y_reverse() +
  theme(panel.grid = element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        panel.border = element_blank())
ggsave('H1_bubble_different_order.pdf', width = 10, height = 10, units = c("cm"))



# figure effector sets
md_ES <- distinct_ES
md_ES$ES_n <- paste(md_ES$ES, ' (',md_ES$n, ')', sep='')
y_order <- (md_ES$ES_n)
md_ES$number_of_eff <- NULL
md_ES$ES <- NULL
md_ES$n <- NULL
data_ES_fig <- melt(data = md_ES, id = 'ES_n')
colnames(data_ES_fig) <- c('ES', 'Locus', 'Effector')
data_ES_fig$ES <- with(data_ES_fig, factor(ES, levels = y_order))

ggplot(data_ES_fig, aes(x = factor(Locus, level = loci), forcats::fct_rev(ES), fill = Effector)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(breaks = effectors, values = mycols) +
  theme(legend.position = 'none') +
  xlab('Effector locus') +
  ylab('Effector set') +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1, vjust= 1)) +
  theme(axis.title.x = element_text(size = 10, angle = 0, hjust=0.5, vjust= -1)) +
  theme(axis.text.y = element_text(size = 10, angle = 0, hjust=1, vjust=0.5)) +
  theme(axis.title.y = element_text(size = 10, angle = 90, hjust= 0.5, vjust= 1)) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.ticks.x = element_line(size = 0.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())
h <- length(md_ES$ES_n)/2
ggsave('H1_ES.pdf', width = 7, height = h, units = c("cm"))




# ------------- H2 ---------------
md_H1 <- subset(md, select = H2)
#### distinct effector sets
distinct_ES <- md_H1

# identify distinct ES and count how often they occur
distinct_ES <- distinct_ES %>% group_by(across()) %>% count

# give name according to prevalence
distinct_ES <- distinct_ES[order(desc(distinct_ES$n)),]
distinct_ES$ES <- paste('H2_ES', as.character(1:nrow(distinct_ES)), sep = '')
distinct_ES_H2 <- distinct_ES

write.csv(distinct_ES, 'H2_distinct_ES.csv', row.names = F)


# figure: bubble chart
# info: https://cran.r-project.org/web/packages/packcircles/vignettes/progressive_packing.html
md_bubble <- data.frame(distinct_ES$ES, distinct_ES$n)
colnames(md_bubble) <- c('ES', 'n')

bubble_v <- c(md_bubble$n)
packing <- circleProgressiveLayout(bubble_v)
md_bubble_fig <- circleLayoutVertices(packing, npoints = length(bubble_v))

ggplot(md_bubble_fig) +
  geom_polygon(aes(x,y,group = id)) +
  coord_equal() +
  theme_bw() + 
  scale_y_reverse() +
  theme(panel.grid = element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        panel.border = element_blank())
ggsave('H2_bubbles.pdf',  width = 10, height = 10, units = c("cm"))

# to get different layouts
ncircles <- length(bubble_v) 
packings <- lapply(1,
                   function(i) { 
                     x <- sample(bubble_v, ncircles) 
                     circleProgressiveLayout(x) 
                   })

packings <- do.call(rbind, packings)
npts <- length(bubble_v)
dat_gg <- circleLayoutVertices(packings, npoints = npts) 

ggplot(data = dat_gg, aes(x, y)) + 
  geom_polygon(aes(group = id),
               fill = 'maroon', color = 'white') +
  coord_equal() +
  theme_bw() + 
  scale_y_reverse() +
  theme(panel.grid = element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        panel.border = element_blank())
ggsave('H2_bubble_different_order.pdf', width = 10, height = 10, units = c("cm"))



# figure effector sets
md_ES <- distinct_ES
md_ES$ES_n <- paste(md_ES$ES, ' (',md_ES$n, ')', sep='')
y_order <- (md_ES$ES_n)
md_ES$number_of_eff <- NULL
md_ES$ES <- NULL
md_ES$n <- NULL
data_ES_fig <- melt(data = md_ES, id = 'ES_n')
colnames(data_ES_fig) <- c('ES', 'Locus', 'Effector')
data_ES_fig$ES <- with(data_ES_fig, factor(ES, levels = y_order))

ggplot(data_ES_fig, aes(x = factor(Locus, level = loci), forcats::fct_rev(ES), fill = Effector)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(breaks = effectors, values = mycols) +
  theme(legend.position = 'none') +
  xlab('Effector locus') +
  ylab('Effector set') +
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust=1, vjust= 1)) +
  theme(axis.title.x = element_text(size = 10, angle = 0, hjust=0.5, vjust= -1)) +
  theme(axis.text.y = element_text(size = 10, angle = 0, hjust=1, vjust=0.5)) +
  theme(axis.title.y = element_text(size = 10, angle = 90, hjust= 0.5, vjust= 1)) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.ticks.x = element_line(size = 0.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())
h <- length(md_ES$ES_n)/2
ggsave('H2_ES.pdf', width = 12, height = h, units = c("cm"))





# ------------- H3 ---------------
md_H1 <- subset(md, select = H3)
#### distinct effector sets
distinct_ES <- md_H1

# identify distinct ES and count how often they occur
distinct_ES <- distinct_ES %>% group_by(across()) %>% count

# give name according to prevalence
distinct_ES <- distinct_ES[order(desc(distinct_ES$n)),]
distinct_ES$ES <- paste('H3_ES', as.character(1:nrow(distinct_ES)), sep = '')
distinct_ES_H3 <- distinct_ES
write.csv(distinct_ES, 'H3_distinct_ES.csv', row.names = F)


# figure: bubble chart
# info: https://cran.r-project.org/web/packages/packcircles/vignettes/progressive_packing.html
md_bubble <- data.frame(distinct_ES$ES, distinct_ES$n)
colnames(md_bubble) <- c('ES', 'n')

bubble_v <- c(md_bubble$n)
packing <- circleProgressiveLayout(bubble_v)
md_bubble_fig <- circleLayoutVertices(packing, npoints = length(bubble_v))

ggplot(md_bubble_fig) +
  geom_polygon(aes(x,y,group = id)) +
  coord_equal() +
  theme_bw() + 
  scale_y_reverse() +
  theme(panel.grid = element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        panel.border = element_blank())
ggsave('H3_bubbles.pdf',  width = 10, height = 10, units = c("cm"))

# to get different layouts
ncircles <- length(bubble_v) 
packings <- lapply(1,
                   function(i) { 
                     x <- sample(bubble_v, ncircles) 
                     circleProgressiveLayout(x) 
                   })

packings <- do.call(rbind, packings)
npts <- length(bubble_v)
dat_gg <- circleLayoutVertices(packings, npoints = npts) 

ggplot(data = dat_gg, aes(x, y)) + 
  geom_polygon(aes(group = id),
               fill = 'maroon', color = 'white') +
  coord_equal() +
  theme_bw() + 
  scale_y_reverse() +
  theme(panel.grid = element_blank(), 
        axis.text=element_blank(),
        axis.ticks=element_blank(), 
        axis.title=element_blank(),
        panel.border = element_blank())
ggsave('H3_bubble_different_order.pdf', width = 10, height = 10, units = c("cm"))



# figure effector sets
md_ES <- distinct_ES
md_ES$ES_n <- paste(md_ES$ES, ' (',md_ES$n, ')', sep='')
y_order <- (md_ES$ES_n)
md_ES$number_of_eff <- NULL
md_ES$ES <- NULL
md_ES$n <- NULL
data_ES_fig <- melt(data = md_ES, id = 'ES_n')
colnames(data_ES_fig) <- c('ES', 'Locus', 'Effector')
data_ES_fig$ES <- with(data_ES_fig, factor(ES, levels = y_order))

ggplot(data_ES_fig, aes(x = factor(Locus, level = loci), forcats::fct_rev(ES), fill = Effector)) +
  geom_tile(width = 0.9, height = 0.9) +
  scale_fill_manual(breaks = effectors, values = mycols) +
  theme(legend.position = 'none') +
  xlab('Effector locus') +
  ylab('Effector set') +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust=1, vjust= 1)) +
  theme(axis.title.x = element_text(size = 10, angle = 0, hjust=0.5, vjust= -1)) +
  theme(axis.text.y = element_text(size = 8, angle = 0, hjust=1, vjust=0.5)) +
  theme(axis.title.y = element_text(size = 10, angle = 90, hjust= 0.5, vjust= 1)) +
  theme(axis.ticks.y = element_blank()) +
  theme(axis.ticks.x = element_line(size = 0.25)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank())
ggsave('H3_ES.pdf', width = 5, height = 6, units = c("cm"))

