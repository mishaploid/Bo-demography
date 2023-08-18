library(tidyverse)

results <- read.csv("~/Documents/bo-demography/models/selective_sweeps/raisd_full_results.csv",
                    header = TRUE)

results_corrected <- results %>%
  filter(chr %in% c('C1', 'C5', 'C6', 'C9') |
           chr %in% 'C2' & !window_start %in% 39330640:40574312 |
           chr %in% 'C3' & !window_start %in% 46416848:48392988 |
           chr %in% 'C4' & !window_start %in% 21410556:22776112 |
           chr %in% 'C7' & !window_start %in% c(28739:2504186, 13086101:14378995) |
           chr %in% 'C8' & !window_start %in% 12104131:13969358)

# # plot mu statistic by chromosome
results %>%
  group_by(morph) %>%
  filter(chr == 'C8',
         morph %in% c('italica', 'botrytis')) %>% 
  ungroup() %>% 
  mutate(morph = fct_recode(morph,
                            'broccoli' = 'italica',
                            'cauliflower' = 'botrytis')) %>% 
  mutate(col = (mu > quantile(mu, c(0.999)))) %>%
  ungroup() %>%
  ggplot(., aes(position, mu, col = col)) +
  geom_point(size = 2) +
  scale_color_manual(values = c('black', 'orange')) +
  facet_grid(morph ~ chr) +
  theme_bw() +
  ylab(expression(mu)) + 
  theme(legend.position = 'none',
        text = element_text(size = 24),
        axis.text = element_text(size = 14))

ggsave("~/Desktop/sweep_test.png", height = 5, width = 5)
