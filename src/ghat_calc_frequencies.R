################################################################################
## Prep allele frequency data for Ghat 
## S. Turner-Hissong 
## 9 December 2019
################################################################################

library(Ghat)
library(parallel)

# load data
load("~/Desktop/ghat/ghat_data.RData")

pops <- list() 

for (morph in unique(pheno_corrected$morphotype)) {
  print(morph) 
  samps <- filter(pheno_corrected, morphotype %in% morph) %>%
    select(id) 
  pops[[morph]] <- gt_imp[samps$id,]
}

pops[['hilarionis']] <- NULL

## frequencies

frequencies <- map_df(pops, .f = ~colMeans(.x)/2) 

frequencies %>%
  gather(key = morphotype, value = frequency) %>%
  ggplot(., aes(frequency)) +
  geom_histogram() + 
  facet_wrap(~morphotype) 

# pairwise differences in allele frequencies
# solution from: https://stackoverflow.com/questions/28187402/calculate-pairwise-difference-between-each-pair-of-columns-in-dataframe
combs <- outer(colnames(frequencies),
               colnames(frequencies),
               paste,
               sep = "_")

# indx1 <- which(lower.tri(combs, diag = TRUE))

res <- outer(1:ncol(frequencies),
             1:ncol(frequencies),
             function(x,y) 
               frequencies[,x] - frequencies[,y]) 

colnames(res) <- combs

# remove diagonal
res <- res[, which(colSums(res) != 0)]

# plot difference in allele frequency for each pair
res %>%
  summarize_all(mean) %>% 
  gather("morphotypes", "value") %>%
  separate(morphotypes, into = c("morph1", "morph2")) %>% 
  ggplot(., aes(fct_reorder(morph1, -value), 
                fct_reorder(morph2, -value))) +
  geom_tile(aes(fill = value)) + 
  scale_fill_viridis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle("Difference in allele frequency") +
  xlab("Population 1") +
  ylab("Population 2")

################################################################################
## Use LD decay to determine number of effective markers
################################################################################

map <- data.frame(Name = colnames(gt_imp)) %>%
  separate(Name, '_', into = c("Chromosome", "Position"), remove = FALSE) %>%
  mutate(Chromosome = str_replace(Chromosome, 'C', '')) %>%
  mutate_if(is.character, list(as.integer)) %>%
  mutate_if(is.factor, as.character)

gt_impv2 <- as.data.frame(gt_imp) %>%
  mutate(X = rownames(gt_imp)) %>%
  select(X, everything()) %>%
  mutate_if(is.numeric, as.integer)

ld <- ld_decay2(gen = gt_impv2, map = map,
                max_win_snp = 2000, 
                max.chr = 9,
                cores = 1, 
                max_r2 = 0.03)

################################################################################
## Run Ghat
################################################################################

Ghat_area <- apply(res, 2, FUN = function(x) {
  Ghat(effects = effects_corrected$area,
       change = x,
       method = "scale",
       perms = 1000,
       plot = "Ghat",
       num_eff = 864)
}
)

### try for uncorrected data 

Ghat_area_uncorrected <- apply(res, 2, FUN = function(x) {
  Ghat(effects = effects_uncorrected$area,
       change = x,
       method = "scale",
       perms = 1000,
       plot = "Ghat",
       num_eff = 864)
}
)