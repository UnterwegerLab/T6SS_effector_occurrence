# clear R's brain
rm(list = ls())

# load packages
library(dplyr)        # data manipulation package
library(ggplot2)      # figure making
library(gridExtra)    # multiple grid-based plots can be arranged on a page
library(reshape2)     # transform data between wide and long formats
library(stringr)      # wrappers for common string operations


db <- "allTminus"


#####PATH SETTING
#path to db results folder
pathDb <- ""
#path to RawData from blastn search
pathRaw <- 'raw_data/'
#path with RawData plus headings
pathHead <- 'headers/'
#path fo saving of indivudal files (.csv, .pdf) from R analysis
pathR <- 'results/'



#####DATA LOADING
#patient information
myisolates <- read.csv(file = 'sequences_IDs_allTminus.csv', header=TRUE)

blast_results <- myisolates


#blastn raw data (load, give headings, save as new .csv with headings)
add_headers_to_raw_blastn_files <- function(raw_data_file) {
  raw_data <- read.csv(file = paste(pathRaw, raw_data_file, sep=''), header = FALSE)
  colnames(raw_data) <- c("qseqid", "qlen", "qstart", "qend", "sacc", "slen", "sstart", "send", "qseq", "sseq", "length","nident", "pident", "bitscore", "evalue")
  # extract sequence ID
  Sequence_ID <- gsub('\\___.*', '', raw_data$sacc)
  # move column with sequence ID to the beginning
  raw_data <- cbind(Sequence_ID, raw_data)
  write.csv(raw_data, file = paste(pathHead, raw_data_file, sep=''), row.names = FALSE, na="")
}

file_list <- list.files(pathRaw)
for (file_name in file_list) {
  n <- paste(pathRaw, file_name, sep='')
  if(file.size(n) > 0) {
    add_headers_to_raw_blastn_files(file_name)
  }
}

all_files <- c()
file_list <- list.files(pathHead)


######PROCESSING OF DATA
for (n in file_list) {
  #n <- 'allT_PA0092_DK09.csv'
  
  #generate a as in: locus_query
  a <- gsub('.{4}$', '', n)
  a <- gsub(db, '', a)
  a <- gsub('^_', '', a)
  
  #make locus_list for steps after this for loop
  all_files <- append(all_files, a)
  
  #give output file names
  eocname <- paste(db, "_", a, "_eoc.pdf", sep='')
  figname <- paste(db, "_", a, "_R.pdf", sep='')
  onelocus <- paste(db,"_", a, "_R.csv", sep='')
  
  #load raw data (with headings)
  k <- read.csv(file = paste(pathHead, n, sep=''), header = TRUE)
  
  #calculate flp as in full length pident
  k$flp = (k$pident * (k$length)/(k$qlen + (k$length - (k$qend - k$qstart +1))))
  k$alignment_length_percent = 100 * k$length/(k$qlen + (k$length - (k$qend - k$qstart + 1)))
  
  #calculation of eoc-index and summing up of flp and alignment_length_percent
  x1 <- k %>% group_by(Sequence_ID) %>% summarise(n = n())
  colnames(x1) <- c("Sequence_ID", "hits")
  k <- group_by(k, Sequence_ID)
  k <- dplyr::full_join(k, x1, by = "Sequence_ID")
  k <- k %>% mutate(flp_sum = sum(flp))
  k <- k %>% mutate(alignment_length_percent_sum = sum(alignment_length_percent))
  k <- ungroup(k, Sequence_ID)
  k <- k %>% mutate(eoc=ifelse(k$hits == 1, 0, k$flp/k$alignment_length_percent))
  k <- k %>% mutate(eoc_combined=ifelse(k$hits == 1, 0, k$flp_sum/k$alignment_length_percent_sum))
  
  #make figure with eoc-index and sums (flp_sum, alignment_length_percent_sum)
  k_eoc <- k
  k_eoc <- k_eoc %>% filter(k_eoc$hits > 1)
  ggplot(k_eoc, aes(x=Sequence_ID)) +
    geom_point(aes(y = eoc), size = 4, shape = 16, color = "plum") +
    theme_bw(base_size=20) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
  ggsave(filename=paste(pathR,'/',eocname,sep=""), width = 14, height = 14, units = c("in"))
  
  #write .csv with eoc
  l <- k[(k$hits>1),]
  write.csv(l, file = paste(pathR, db, "_", a, "_only_eoc.csv" , sep=''), row.names = FALSE, na="")
  
  #only best hit per genotype / glycerol no
  k <- group_by(k, Sequence_ID)
  k <- arrange(k, desc(flp))
  k <- top_n (k, 1, flp)
  k <- ungroup(k)
  
  #figure making data.frame
  k <- k %>% mutate(help_column = ifelse(k$hits > 1,
                                         ifelse(k$eoc_combined > 0.95,
                                                100, k$flp),
                                         k$flp))
  if (str_detect(n, "PA5086-88") == TRUE) {k <- k %>% mutate(help_column = ifelse(k$alignment_length_percent == 100, k$help_column, 0))}
  if (str_detect(n, "PA3291") == TRUE) {k <- k %>% mutate(help_column = ifelse(k$alignment_length_percent > 95, k$help_column, 0))}
  if (str_detect(n, "PA3290") == TRUE) {k <- k %>% mutate(help_column = ifelse(k$alignment_length_percent_sum > 95, k$help_column, 0))}
  if (str_detect(a, "STE") == TRUE) {k <- k %>% mutate(help_column = k$flp)}
  
  #make figure: one for each locus
  b <- str_sub(a,-1,-1)
  if (str_detect(b, 'zi') == TRUE) {b <- 'z'}
  if (str_detect(b, 'zz') == TRUE) {b <- 'z'}
  if (str_detect(b, 'za') == TRUE) {b <- 'z'}
  if (str_detect(b, 'zb') == TRUE) {b <- 'z'}
  
  if (str_detect(b, "z") == TRUE) {yesz <- k}
  
  k <- k %>% mutate(flp_figure = k$help_column)
  k$help_column <- NULL
  
  if (str_detect(b, "z") == TRUE) {
    #check if z
    yesz <- yesz %>% mutate(flp_figure = ifelse(yesz$hits == 1,
                                                ifelse(yesz$flp > 90, yesz$help_column, 0),
                                                0))
    yesz$help_column <- NULL}
  
  if(str_detect(b, "z") == TRUE) {write.csv(yesz, file = paste(pathR, onelocus, sep=''), row.names = FALSE, na="")}
  if(str_detect(b, "z") == TRUE) {ggplot(yesz, aes(x=Sequence_ID, y=flp_figure)) +
      geom_point(size=5, shape=15) +
      theme_bw(base_size=20) +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
    ggsave(filename=paste(pathR,'/',figname,sep=""), width = 14, height = 14, units = c("in")) }
  
  if(str_detect(b, "z") == FALSE) {write.csv(k, file = paste(pathR, onelocus, sep=''), row.names = FALSE, na="")}
  if(str_detect(b, "z") == FALSE) {ggplot(k, aes(x=Sequence_ID, y=flp_figure)) +
      geom_point(size=5, shape=15) +
      theme_bw(base_size=20) +
      theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5))
    ggsave(filename=paste(pathR,'/',figname,sep=""), width = 14, height = 14, units = c("in")) }
  
  #merging of results in one data.frame (mydata)
  if(str_detect(b, "z") == FALSE) {m <- data.frame(k$Sequence_ID, k$flp_figure)}
  if(str_detect(b, "z") == TRUE) {m <- data.frame(yesz$Sequence_ID, yesz$flp_figure)}
  colnames(m) <- c("Sequence_ID", a)
  blast_results <- dplyr::full_join(blast_results, m, by = "Sequence_ID")
  
}

write.csv(all_files, file = 'all_files.csv', row.names=FALSE, na='')
write.csv(blast_results, file = paste(db, "_flp.csv", sep=""), row.names = FALSE, na="")



##########################################################
locus_list <- all_files
locus_list <- gsub("\\_.*", "", locus_list)
locus_list <- unique(locus_list)

allkinds <- myisolates

for (n in locus_list) {
  #n <- 'PA14-33970'
  
  #only load data from one locus (from blast_results data.frame)
  m <- paste(n, "_",sep="")
  if(str_detect(n, "PA5086-88") == TRUE) {m <- "PA5086"}
  k <- blast_results %>% dplyr::select(Sequence_ID, starts_with(m))
  k$help_column <- 50
  k$kind <- 1
  k[is.na(k)] <- 1
  
  #specify kind for each isolate
  count <- ncol(k)
  a <- (count - 3 + 1)
  number_of_queries <- seq(2, a)
  column = 0
  k$count_over_80 = 0
  for (o in number_of_queries) {
    column <- count - o
    p <- colnames(k)[column]
    p <- gsub('^.*_', '', p)
    right <- column + 1
    k <- k %>% mutate(kind = ifelse(k[,column] > k$help_column, p, k$kind))
    k <- k %>% mutate(help_column = ifelse(k[,column] > k$help_column, k[,column], k$help_column))
    
    
    k <- k %>% mutate(count_over_80 = ifelse(k[,column] > 80, k$count_over_80 + 1, k$count_over_80))
    
  }
  
  k$help_column <- NULL
  k <- k %>% mutate(kind = ifelse(k$count_over_80 > 1, 'ybdl_mipeg_T6aog_mhfz', k$kind))
  
  # ybdl: yes but different locus (PA0260, PA0262)
  # mipeg: multiple immunity protein-encoding genes (PA0092, PA3290)
  # T6aog: T6 and other genes (PA0822, PA0821)
  # mhfz: multiple hits for different z (PA3487)
  
  k$count_over_80 <- NULL
  
  # check for ybdl (yes but at different locus) (i.e. T_1221 PA0260)
  
  # modify for PA0822 (2023-01-19)
  if (str_detect(n, "PA0822") == TRUE) {
    k <- k %>% mutate(kind = ifelse(k$kind == 'ybdl_mipeg_T6aog_mhfz',
                               ifelse(k$PA0822_PAO1 == 1, 'STEX', k$kind),
                               k$kind))
  }
  
  # modify for PA14-33970 (2023-01-19)
  if (str_detect(n, "33970") == TRUE) {
    k <- k %>% mutate(kind = ifelse(k$kind == 'ybdl_mipeg_T6aog_mhfz',
                                    ifelse(k$`PA14-33970_DK01` == 1, 'STEX', k$kind),
                                    k$kind))
  }
  
  #write .csv
  write.csv(k, file = paste(pathR, db, "_", n, ".csv", sep=""), row.names = FALSE, na="")
  
  #combine in one big file
  l <- data.frame(k$Sequence_ID, k$kind)
  colnames(l) <- c("Sequence_ID", paste(n))
  allkinds <- dplyr::full_join(allkinds, l, by = "Sequence_ID")
}


write.csv(allkinds, file = paste(db, "_allkinds.csv", sep=""), row.names = FALSE, na="")

