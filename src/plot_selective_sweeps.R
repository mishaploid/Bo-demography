library(tidyverse)
library(GenomicRanges)

foo <- read.csv("~/Documents/bo-demography/models/selective_sweeps/sig_windows.csv")

foo %>%
  # filter(chr == 'C1') %>% 
  ggplot(.) + 
  geom_segment(aes(x = window_start, xend = window_end,
                   y = morph, yend = morph),
               stat = "identity",
               size = 4,
               alpha = 0.25) +
  facet_grid(~chr) + 
  theme_bw()

testdf <- foo

gr = GRanges(seqnames = testdf$chr,
             IRanges(testdf$window_start,
                     testdf$window_end),
             morph = testdf$morph,
             strand = "*") 

gr_split <- reduce(split(gr, elementMetadata(gr)$morph))
gr_df <- as.data.frame(gr_split)

gr_split

mask <- read.table("~/Documents/bo-demography/data/processed/mappability_masks/scratch/combined_mask.bed",
                   header = FALSE, col.names = c("chr", "start", "stop"))

mask_gr <- GRanges(seqnames = mask$chr,
                    IRanges(mask$start,
                            mask$stop),
                    strand = "*") 

mask2 <- mask_gr[width(mask_gr) < 40000000] 

mask2 <- as.data.frame(reduce(mask2))

gr_df %>%
  # filter(seqnames %in% c('C1', 'C5', 'C6', 'C9') |
  #        seqnames %in% 'C2' & !start %in% 3933064:40574312 |
  #        seqnames %in% 'C3' & !start %in% 46416848:48392988 |
  #        seqnames %in% 'C4' & !start %in% 21410556:22776112 |
  #        seqnames %in% 'C7' & !start %in% c(28739:2504186, 13086101:14378995) |
  #        seqnames %in% 'C8' & !start %in% 12104131:13969358) %>%
  mutate(group_name = fct_relevel(group_name,
                             'alboglabra', 'italica', 'botrytis', 'gongylodes', 'capitata', 'sabauda', 'gemmifera',
                             'costata', 'medullosa', 'palmifolia', 'ramosa', 'sabellica', 'viridis',
                             'oleracea', 'cretica', 'incana', 'insularis', 'rupestris', 'villosa'),
         type = ifelse(group_name %in% c('cretica', 'oleracea', 'incana',
                                         'insularis', 'rupestris', 'villosa'), 
                       paste('wild'), paste('domesticated'))) %>%
  ggplot(., aes(y = group_name, yend = group_name,
                  x = start, xend = end,
                col = type)) +
  scale_colour_manual(values = c("#009E73",
                                 "#E69F00")) + 
  geom_segment(size = 1.5) +
  geom_rect(data = mask2, aes(ymin = -Inf, ymax = Inf, xmin = start, xmax = end),
            inherit.aes = FALSE, fill = 'black') +
  facet_grid(seqnames ~ ., space = "free", scale = "free") +
  theme(legend.position = "none") 

ggsave("~/Documents/bo-demography/models/selective_sweeps/figures/significant_windows.png",
       height = 15, width = 15)


head(gr_df)

foo <- countOverlaps(gr_split[seqnames(gr_split) == 'C1'], 
              type = "any", 
              drop.self = FALSE)



# test <- fuzzy_inner_join(foo, foobar, 
#                          by=c("V1"="WindowStart", "V1"="WindowStop"),
#                          match_fun=list(`>=`, `<=`))



# library(IRanges)

# test <- IRanges(gem$WindowStart, gem$WindowStop)
# test2 <- IRanges(bot$WindowStart, bot$WindowStop)
# 
# c1 <- sig_win %>% 
#   group_by(chr) %>% 
#   arrange(WindowStart, WindowStop) 

sig_win <- as.data.frame(gr_split)

overlaps <- list() 

for (i in 1:9) {
  tmp <- sig_win %>%
    filter(seqnames %in% paste0('C', i)) %>%
    arrange(start, end)
  ranges <- IRanges(tmp$start, tmp$end)
  overlaps[[i]] <- data.frame(morph = tmp$group_name,
                              start = tmp$start,
                              end = tmp$end,
                              chr = paste0('C', i), 
                              overlaps = countOverlaps(ranges))
}

library(data.table)

overlap_df <- rbindlist(overlaps)

library(viridisLite)

overlap_df2 <- overlap_df %>% 
  mutate(overlap = ifelse(overlaps >= 1, 1, 0)) %>%
  select(-overlaps) %>% 
  # rownames_to_column() %>%
  spread(key = morph, value = overlap) 

seqs <- GRanges(seqnames = overlap_df$chr,
             IRanges(overlap_df$start,
                     overlap_df$end),
             morph = overlap_df$morph,
             strand = "*") 

sig_regions <- as.data.frame(reduce(seqs))

foo <- overlap_df2 %>% 
  group_by(chr) %>% 
  arrange(chr, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) > 
                              cummax(as.numeric(end)))[-n()])) %>% 
  gather(key = morph, value = overlap, -start, -end, -chr, -indx) %>%
  # filter(chr == 'C9', indx == 10) %>%
  group_by(chr, indx) %>% 
  na.omit() %>% 
  filter(!duplicated(morph)) %>%
  mutate(count = n()) %>%
  ungroup() %>%
  group_by(chr) %>%
  filter(!duplicated(indx))
  
foo %>% 
  mutate(morph = fct_relevel(morph,
                                  'alboglabra', 'italica', 'botrytis', 'gongylodes', 'capitata', 'sabauda', 'gemmifera',
                                  'costata', 'medullosa', 'palmifolia', 'ramosa', 'sabellica', 'viridis',
                                  'oleracea', 'cretica', 'incana', 'insularis', 'rupestris', 'villosa'),
         type = ifelse(morph %in% c('cretica', 'oleracea', 'incana',
                                         'insularis', 'rupestris', 'villosa'), 
                       paste('wild'), paste('domesticated'))) %>% 
  ggplot(aes(count, fill = cut(count, 18))) +
  geom_histogram(bins = 19) +
  scale_fill_viridis_d(option = "magma") + 
  theme_classic() + 
  theme(text = element_text(size = 24),
        legend.position = "none",
        panel.background = element_rect(fill='gray90')) +
  xlab('Number of morphotypes') 
  
  
ggsave("~/Documents/bo-demography/models/selective_sweeps/shared_sweep_histogram.png",
       height = 5, width = 5)

chr_length <- mask %>% 
  group_by(chr) %>%
  summarize_all(.funs = c(min, max))

library(viridisLite)

ggplot() +
  # segments for full chromosome length
  geom_segment(data = chr_length,
               aes(x = start_fn1, 
                   xend = stop_fn2,
                   y = fct_rev(chr), 
                   yend = fct_rev(chr)),
               color = "gray90",
               size = 14) +
  # segments for selective sweeps
  geom_segment(data = foo,
               aes(x = start, 
                   xend = end, 
                   y = fct_rev(chr), 
                   yend = fct_rev(chr),
                   color = count),
               size = 14) +
  scale_color_viridis_c(option = "magma",
                        'morphotypes') + 
  theme_classic() +
  theme(text = element_text(size = 22),
        axis.text.x = element_text(size = 14),
        legend.position = "top") +
  xlab('Position') + 
  ylab('Chromosome')

ggsave("~/Documents/bo-demography/models/selective_sweeps/figures/genomewide_sweeps_by_morph.png",
       height = 6, width = 9)

ggplot(overlap_df, aes(overlaps, fill = morph)) +
  geom_histogram() +
  scale_fill_viridis_d('Morphotype') + 
  # facet_wrap(~chr) +
  xlab('Number of overlapping regions') +
  theme_classic() +
  theme(text = element_text(size = 24))

ggsave('~/Documents/bo-demography/models/selective_sweeps/figures/overlapping_windows.png', height = 6, width = 7)

length(overlap_df$WindowStart)

library(IntervalSurgeon)

join(x, y, output = 'indices')

sig_win %>%
  filter(morph %in% c('alboglabra', 'italica', 'botrytis', 'gongylodes', 
                      'capitata', 'sabauda', 'gemmifera', 'costata', 'medullosa',
                      'palmifolia', 'ramosa', 'sabellica', 'viridis')) %>% 
ggplot(., aes(x = WindowStart, xend = WindowStop, y = morph, yend = morph)) +
  geom_segment(size = 2) + 
  facet_grid(chr ~ .) +
  theme_minimal()

ggsave("~/Documents/bo-demography/models/selective_sweeps/figures/overlapping_regions.png", height = 12, width = 20)


foo <- windows %>%
  filter(chr %in% 'C1') %>% 
  group_by(morph) %>%
  do(Intervals(c(WindowStart, WindowStop)))
  unite(temp, WindowStart, WindowStop) %>% 
  spread(key = morph, value = temp)
  
  

library(intervals)
test <- Intervals(c(gem$WindowStart, gem$WindowStop))
test2 <- Intervals(c(bot$WindowStart, bot$WindowStop))
interval_overlap(test, test2)

gem <- gr[gr$morph == 'gemmifera']
bot <- gr[gr$morph == 'botrytis']

multiOverlap(gem, bot)

gr_test <- gr[seqnames(gr) == 'C1']
gr_test

test <- countOverlaps(gr_test)
test2 <- countOverlaps(gr_test, gr_test[test==1])


foverlaps(gem, bot, type="within", nomatch=0L, which=TRUE)[, .N, by=c(WindowStart, WindowStop)]
