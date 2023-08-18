library(tidyverse)
library(cowplot)

cols <- c(colorRampPalette(c("black", "gray"))(4),
          colorRampPalette(c("#004949", "#66c2a4", "#ccece6"))(6),
          "#490092", "#AA4499", "#DDCC77", "#006ddb", "#882255", "#924900") 

names(cols) <- c("Brassica_cretica", "Brassica_incana", "Brassica_rupestris", "Brassica_villosa",
                 "tronchuda cabbage", "marrow cabbage", "lacinato kale", "perpetual kale", "collards", "curly kale",
                 "broccoli", "cauliflower", "cabbage", "kohlrabi", "Chinese kale", "Brussels sprouts")

new_names <- c("Brassica oleracea" = "oleracea",
               "tronchuda cabbage" = "costata",
               "marrow cabbage" = "medullosa",
               "perpetual kale" = "ramosa",
               "lacinato kale" = "palmifolia",
               "curly kale" = "sabellica",
               "collards" = "viridis",
               "Chinese kale" = "alboglabra",
               "cabbage" = "capitata",
               "Brussels sprouts" = "gemmifera",
               "kohlrabi" = "gongylodes",
               "broccoli" = "italica",
               "cauliflower" = "botrytis")

smc_results <- read_csv("reports/smc_cv_no_timepoints_results.csv") %>% 
  mutate(population = ifelse(grepl("Brassica_|oleracea", label),
                             "wild/weedy",
                             ifelse(grepl("medullosa|palmifolia|ramosa|sabellica|viridis|costata|alboglabra", label),
                             "kales",
                             "cultivated")),
         label = fct_recode(label, !!!new_names),
         label = fct_relevel(label,
                             "curly kale", after = Inf))

# split data for each "facet"
smc_results <- split(smc_results, f = smc_results$population)

# plot for the first "facet"
p1 <- smc_results$`wild/weedy` %>% 
  # filter(population %in% "wild/weedy") %>% 
  ggplot(., aes(x = x, 
                y = y,
                color = label)) + 
  geom_line(size = 1.5,
            alpha = 0.8) +
  scale_color_manual(name = "",
                     values = cols) + 
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 4),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(400, 10^5.5)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x, n = 2),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                limits = c(10^3.2, 10^5.5)) +
  facet_wrap(population ~ ., ncol = 1) +
  annotation_logticks() +
  xlab("Generations ago") +
  ylab("Effective population size") +
  theme_classic() + 
  theme(legend.justification = "left") 

p2 <- p1 %+% smc_results$kales 

p3 <- p1 %+% smc_results$cultivated

plot_grid(p1, p2, p3,
          ncol = 1,
          align = "v",
          labels = "AUTO")

ggsave("reports/figures/Figure2.png",
       height = 8, width = 6)
