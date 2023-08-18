library(tidyverse)
library(vroom)
library(RColorBrewer)

# specify color palette (replace number in parentheses with max k value)
cols <- colorRampPalette(brewer.pal(8, "Dark2"))(13)

# load sample ids ---------------------------------------------------------
# Want to make sure samples are in correct order & combine with population ids
# I have population ids in a separate file, prob should just replace the fam file

# fam file with sample order
# this is a .fam formatted file from PLINK 
fam <- read.table("biallelic_snps.geno10.maf05.ldpruned.fam") %>%
  select(V2) %>%
  rename(sample = V2)

# population ids for each sample 
# text file with two columns (sample and population)
groups <- read.table("../sample_ids.txt") %>%
  rename(sample = V1,
         population = V2) %>% 
  left_join(fam, ., 
            by = "sample") 

# load admixture results --------------------------------------------------

# list admixture output files 
# assumes files are located in current working directory
afiles <- list.files(path = ".", 
                     pattern = "*.Q$",
                     full.names = TRUE) 

# read in admixture results to a single df 
q_df <- map_df(afiles, ~vroom(.x,
                              col_names = FALSE,
                              id = "k")) %>%
  # may need to adjust this - wanted to extract k value from file name
  mutate(k = word(k, 6, sep = "[.]")) %>%
  # combine results with population ids 
  cbind(groups, .) %>% 
  # convert to long format 
  gather(key = cluster,
         value = prop,
         -k, -sample, -population) %>%
  # remove missing values for lower values of k
  drop_na() 


# sort samples by proportion ----------------------------------------------
# https://stackoverflow.com/questions/41679888/how-to-group-order-data-in-r-for-a-barplot

q_ordered <- q_df %>%
  group_by(sample) %>%
  # determine most likely assignment 
  mutate(likely_assignment = population[which.max(prop)],
         assignment_prob = max(prop)) %>% 
  # sort samples by assignment 
  arrange(likely_assignment, desc(assignment_prob)) %>%
  ungroup() %>%
  # make sure samples are ordered
  mutate(sample = fct_inorder(sample)) 


# plot the results! -------------------------------------------------------

q_ordered %>% 
<<<<<<< HEAD
  ggplot(., aes(sample, propx, fill = cluster)) +
=======
  ggplot(., aes(sample, prop, fill = cluster)) +
>>>>>>> 86e618b2156c55942bad0b79fb75df8e578cc3a0
  # stacked barplot
  geom_bar(stat = "identity") +
  # use custom colors
  scale_fill_manual(values = cols) + 
  # rotate x-axis labels 
  scale_x_discrete(guide = guide_axis(angle = 90)) + 
  # facet by k value and population id
  facet_grid(k ~ population,
             scales = "free_x",
             space = "free") +
  theme_bw() +
  # rotate strip text for easier reading 
  theme(strip.text.x = element_text(angle = 90,
                                    hjust = 0),
        axis.text.x = element_text(size = 4),
        # reduce spacing between panels 
        panel.spacing.x = unit(0.001, "lines"),
        legend.position = "none")

# save the output 
ggsave("admixture_results.pdf",
       height = 18,
       width = 16)