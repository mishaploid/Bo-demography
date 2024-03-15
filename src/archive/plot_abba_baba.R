### ABBA-BABA plotting 
library(tidyverse)

results <- read.csv("~/Documents/bo-demography/models/ABBA_BABA/results/sabellica_palmifolia.csv",
                    header = TRUE) %>%
  mutate(W = fct_recode(W,
                        `Chinese kale` = 'alboglabra',
                        `savoy cabbage (wirsing)` = 'sabauda',
                        kohlrabi = 'gongylodes',
                        collards = 'viridis',
                        cabbage = 'capitata',
                        `tronchuda kale` = 'costata',
                        broccoli = 'italica',
                        `perpetual kale` = 'ramosa',
                        `Brussels sprouts` = 'gemmifera',
                        `marrow cabbage` = 'medullosa',
                        cauliflower = 'botrytis'))


ggplot(results, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr),
                size = 1.5) +
  theme_minimal() + 
  theme(# axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        text = element_text(size = 24)) +
  xlab("H1") +
  coord_flip()

ggsave("~/Desktop/ghat/abba_baba.png", height = 4.5, width = 10)
