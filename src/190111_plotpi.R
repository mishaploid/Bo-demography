library(tidyverse)

diversity <- data_frame(filename = list.files(path = "~/Documents/bo-demography/reports/neutrality", pattern = "*pestPG", full.names = TRUE)) %>%
  mutate(file_contents = map(filename, ~ read_table2(file.path(.), col_names = TRUE))) %>%
  unnest() %>%
  mutate(filename = gsub("(.*)/", "", filename),
         pop = gsub(".pestPG", "", filename),
         pi = tP/nSites,
         pop = factor(pop, levels = c("wild", "sabellica", "acephala", "alboglabra",  "gongylodes", "capitata", "italica", "botrytis", "gemmifera")))

diversity <- diversity %>%
  mutate(pop = recode_factor(pop, sabellica = "curly kale",
                             acephala = "ornamental kale",
                             alboglabra = "Chinese kale",
                             gongylodes = "kohlrabi",
                             capitata = "cabbage",
                             italica = "broccoli",
                             botrytis = "cauliflower",
                             gemmifera = "Brussels sprouts"),
         pop = factor(pop, levels = c("wild", "ornamental kale", "curly kale",
                                      "Chinese kale", "kohlrabi", "cabbage", "broccoli",
                                      "cauliflower", "Brussels sprouts")))

# set the color palette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "black")
colorset <- c("wild" = "black", "cabbage" = "#FFD700", "kohlrabi" = "#56B4E9", 
              "cauliflower" =  "#E69F00", "broccoli" = "#009E73", "Chinese kale" = "#D55E00", 
              "Brussels sprouts" = "#0072B2", "ornamental kale" = "#999999", "curly kale" = "#CC79A7")

ggplot(diversity, aes(reorder(pop, -pi), pi, fill = pop)) +
  geom_boxplot(size = 1) +
  scale_fill_manual(values = colorset) +
  theme_bw(base_size = 24) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(),
        legend.position = "none") +
  ylab(expression(pi))

ggsave("~/Documents/bo-demography/pi.png", height = 5, width = 6)

ggplot(diversity, aes(reorder(pop, -pi), Tajima, fill = pop)) +
  geom_boxplot(size = 1) +
  scale_fill_manual(values = colorset) +
  theme_bw(base_size = 24) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(),
        legend.position = "none") +
  ylab("Tajima's D")

ggsave("~/Documents/bo-demography/TajimasD.png", height = 5, width = 6)

mean_pi <- diversity %>% group_by(pop) %>% summarize_all(mean)

tajd <- data_frame(filename = list.files(path = "../../Bo-popgen/neutrality_stats/TajimasD", pattern = "*Tajima.D", full.names = TRUE)) %>%
  mutate(file_contents = map(filename, ~ read_table2(file.path(.), col_names = TRUE))) %>%
  unnest() %>%
  mutate(filename = gsub("(.*)/", "", filename),
         pop = gsub(".Tajima.D", "", filename)) %>%
  left_join(., mean_pi, by = "pop")

ggplot(tajd, aes(reorder(pop, -PI), TajimaD, fill = pop)) +
  geom_boxplot(size = 1) +
  scale_fill_manual(values = colorset) +
  theme_bw(base_size = 24) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(),
        legend.position = "none") +
  ylab("Tajima's D")

ggsave("~/Documents/bo-demography/TajD.png", height = 6, width = 8)

