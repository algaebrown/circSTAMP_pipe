library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

# Read in data
counts = read_tsv(args[1])
replicates = unlist(strsplit(args[2], ','))
outf = args[3]

# compute rho
get_variance <- function(sample_pairs){
   rep1 = paste0('BSJ-',sample_pairs[1])
   rep2 = paste0('BSJ-',sample_pairs[2])
    
    
               
    print(paste('estimating variance between: ',rep1, rep2))
    
    counts$total_reads=counts[[rep1]]+counts[[rep2]]
    sub_counts = counts[counts[[rep1]]>0 & counts[[rep2]]>0, ]
    
    betabinom_fit = VGAM::vglm(cbind(sub_counts[[rep1]], sub_counts[[rep2]]) ~ 1 ,
                           VGAM::betabinomial(), trace = TRUE) # rho as a function to total reads
    print(VGAM::coef(betabinom_fit, matrix = TRUE))
    betabinom_coef = betabinom_fit %>% (VGAM::coef) %>% as_tibble(rownames="coef") %>% transmute(coef = c("mu","rho"), value)%>% pivot_wider(names_from=coef,values_from=value) %>%mutate(rep1 = rep1, rep2 = rep2)
    
    return(betabinom_coef)
} 

coef = lapply(combn(replicates, 2, simplify=FALSE), get_variance) %>% bind_rows
write_tsv(coef, outf)