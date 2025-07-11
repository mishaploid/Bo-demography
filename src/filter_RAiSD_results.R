#########################################
## Combine and filter RAiSD results
#########################################

# This script reads in results from RAiSD for multiple subspecies, chromosmoes,
# and window sizes. The results are filtered by group to select the top 1% 
# of Mu values. The filtered results are then exported to a tab-delimited file
# for each combination of subspecies and window size. 

library(tidyverse)

# list RAiSD output files and read into a data frame 
sweeps <- list.files('models/raisd',
                     pattern = '^RAiSD_Report.*',
                     full.names = TRUE) %>% 
  set_names() %>% 
  map_dfr(., 
          ~read_delim(.x, delim = '\t',
                      skip = 1,
                      col_names = c('Position', 'WindowStart', 'WindowEnd', 'FactorsVAR', 'SFS', 'LD', 'Mu')),
          .id = 'filename') 

# filter RAiSD output to include only rows in the top 1% of Mu values for each
# combination of subspecies/chr/window size
sig_sweeps <- sweeps %>% 
  group_by(filename) %>% 
  mutate(q99 = quantile(FactorsVAR, probs = c(0.99))) %>% 
  filter(FactorsVAR >= q99) %>% 
  mutate(filename = gsub('Brassica_', 'Brassica', filename),
         subspecies = word(filename, 5, 5, sep = '/|_|\\.'),
         subspecies = gsub('Brassica', 'Brassica_', subspecies),
         chr = word(filename, 6, 6, sep = '/|_|\\.'),
         window_size = word(filename, 7, 7, sep = '/|_|\\.')) %>% 
  ungroup() %>% 
  select(subspecies, chr, window_size, Position, WindowStart, WindowEnd, FactorsVAR, SFS, LD, Mu, q99) 

# export filtered RAiSD results to individual tab-delimited files 
sig_sweeps %>% 
  group_by(subspecies, window_size) %>% 
  group_walk(~write_delim(.x, paste0('models/raisd/q99_sweeps/', .y$subspecies, '_', .y$window_size, '_q99_sweeps.txt'),
                          delim = '\t'),
             .keep = TRUE)
