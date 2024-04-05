# Quick script to compile and plot results from RAiSD
# S. Turner-Hissong
# 2 January 2020

library(tidyverse)
library(GenWin)

# list files containing results from RAiSD
filelist <- list.files(pattern = 'correctedC',
                       recursive = TRUE)

# read in all files at once using purrr::map
results <- tibble(filename = filelist) %>%
  mutate(contents = map(filename,
                        ~read.table(.,
                                    skip = 1,
                                    header = FALSE,
                                    col.names = c('position',
                                                  'window_start',
                                                  'window_end',
                                                  'mu_VAR',
                                                  'mu_SFS',
                                                  'mu_LD',
                                                  'mu')))) %>%
  unnest() %>%
  separate(filename,
           into = c('morph', 'chr'),
           sep = '/') %>%
  mutate(chr = str_replace(chr, 'corrected', ''),
         morph = fct_relevel(morph,
           'alboglabra', 'italica', 'botrytis', 'gongylodes', 'capitata', 'sabauda', 'gemmifera',
           'costata', 'medullosa', 'palmifolia', 'ramosa', 'sabellica', 'viridis',
           'oleracea', 'cretica', 'incana', 'insularis', 'rupestris', 'villosa'))

# windows <- results %>%
#   group_by(morph, chr) %>%
#   filter(!duplicated(position)) %>%
#   do(splineAnalyze(Y = .$mu,
#              map = .$position,
#              smoothness = 100,
#              plotRaw = TRUE,
#              plotWindows = TRUE,
#              method = 4)$windowData)

results %>%
  filter(chr %in% c('C1', 'C5', 'C6', 'C9') |
           chr %in% 'C2' & !window_start %in% 39330640:40574312 |
           chr %in% 'C3' & !window_start %in% 46416848:48392988 |
           chr %in% 'C4' & !window_start %in% 21410556:22776112 |
           chr %in% 'C7' & !window_start %in% c(28739:2504186, 13086101:14378995) |
           chr %in% 'C8' & !window_start %in% 12104131:13969358) %>%
 group_by(morph) %>%
 mutate(col = (mu > quantile(mu, c(0.999)))) %>%
 ungroup() %>%
 filter(col == TRUE) %>%
 write.csv(., "sig_windows.csv", row.names = FALSE)

# # plot mu statistic by chromosome
results %>%
  filter(chr %in% c('C1', 'C5', 'C6', 'C9') |
           chr %in% 'C2' & !window_start %in% 3933064:40574312 |
           chr %in% 'C3' & !window_start %in% 46416848:48392988 |
           chr %in% 'C4' & !window_start %in% 21410556:22776112 |
           chr %in% 'C7' & !window_start %in% c(28739:2504186, 13086101:14378995) |
           chr %in% 'C8' & !window_start %in% 12104131:13969358) %>%
  group_by(morph) %>%
  mutate(col = (mu > quantile(mu, c(0.999)))) %>%
  ungroup() %>%
  ggplot(., aes(position, mu, col = col)) +
  geom_point(size = 3) +
  scale_color_manual(values = c('black', 'orange')) +
  facet_grid(morph ~ chr,
             scale = 'free',
             space = 'free') +
  theme_bw() +
  theme(legend.position = 'none')

ggsave("figures/mu.png", height = 20, width = 20)
#
# # plot mu_var statistic
#
# results %>%
#   ggplot(., aes(position, mu_VAR)) +
#   geom_point(size = 0.5) +
#   facet_grid(morph ~ chr,
#              scale = 'free',
#              space = 'free') +
#   theme(legend.position = 'none') +
#   theme_bw()
#
# ggsave("figures/mu_VAR.png", height = 20, width = 20)
#
# # plot mu_sfs
#
# results %>%
#   ggplot(., aes(position, mu_SFS)) +
#   geom_point(size = 0.5) +
#   facet_grid(morph ~ chr,
#              scale = 'free_x',
#              space = 'free') +
#   theme(legend.position = 'none') +
#   theme_bw()
#
# ggsave("figures/mu_SFS.png", height = 20, width = 20)
#
# # plot mu_ld
#
# results %>%
#   ggplot(., aes(position, mu_LD)) +
#   geom_point(size = 0.5) +
#   facet_grid(morph ~ chr,
#              scale = 'free_x',
#              space = 'free') +
#   theme(legend.position = 'none') +
#   theme_bw()
#
# ggsave("figures/mu_LD.png", height = 20, width = 20)
#
# results %>%
# group_by(morph) %>%
# mutate(col = (mu > quantile(mu, c(0.999)))) %>%
# ungroup() %>%
# filter(col == TRUE) %>%
# ggplot(., aes(pos, fct_rev(morph), col = col)) +
# geom_point() +
# scale_color_manual(values = c("black", "black")) +
# facet_grid(~chr, scale = "free", space = "free") +
# theme(legend.position = "none")
#
# results2 <- results %>% group_by(morph) %>% mutate(grouped_id = row_number()) %>% spread(morph, mu)
#
# test <- results %>% group_by(chr) %>% mutate(bin = cut(position, breaks = c(-Inf, mean(position), Inf)))
#
# results2 <- test %>% ungroup() %>% group_by(morph) %>% select(morph, chr, bin, sig) %>% mutate(group_id = row_number()) %>% spread(morph, sig)
#
#  bar <- foo %>% ungroup() %>% mutate(count = rowSums(.[,3:21]!=0)) %>% gather(morph, num_sites, -chr, -bin, -count)
