# plot popvae results 

library(data.table)

pd <- fread("../popvae/geno10_maf5_ldpruned_silscores_latent_coords.txt") %>%
  mutate(sampleID = str_replace(sampleID, "0_", ""))

names(pd)[1:2] <- c("LD1", "LD2")

sample_ids <- read.table("../sample_ids.txt",
                         header = FALSE,
                         col.names = c("sample_id",
                                       "population"))

pd <- left_join(pd, sample_ids,
                by = c("sampleID" = "sample_id"))

missingness <- read.table("../../scratch/out.imiss",
                          header = TRUE) %>%
  mutate(sampleID = str_replace(INDV, "0_", ""))

sample_info <- read_csv("../../data/external/Cheng2016_sample_info.csv") %>% 
  mutate(sampleID = paste0("SamC_", str_pad(Index, 3, pad = "0")))


grin_info <- read_csv("~/Desktop/oleracea_location_info/germplasm_data/processed/grin_sample_info.csv") %>% 
  mutate(id = paste0("Bo", str_pad(id, 3, pad = "0"), rep))
ipk_info <- read_csv("~/Desktop/oleracea_location_info/germplasm_data/processed/ipk_sample_info.csv") %>% 
  mutate(id = paste0("Bo", str_pad(id, 3, pad = "0") , rep))

genebank_ids <- c(grin_info$id, ipk_info$id)

devtools::install_github("AndiKur4/MaizePal")

library(MaizePal)

cols <- colorRampPalette(maize_pal("HopiBlue"))(8)
cols <- c(cols, "gray40")

pd %>%
  mutate(id = str_extract(sampleID, "(Bo)(.*?)(?=_)")) %>% 
  mutate(data_source = ifelse(str_detect(sampleID, "Bo"),
                         "Mabry",
                         ifelse(str_detect(sampleID, "Sam"),
                                "Cheng",
                                "Kioukis")),
         population = ifelse(population %in% "sabauda",
                             "capitata", 
                             paste0(population)),
         population = ifelse(grepl("Brassica_", population),
                             "wild", 
                             paste0(population)),
         population = ifelse(grepl("medullosa|palmifolia|ramosa|sabellica|viridis|costata", population),
                                   "kale",
                                   paste0(population)),
         population = fct_recode(population,
                                 "Chinese kale" = "alboglabra",
                                 "cauliflower" = "botrytis",
                                 "cabbage" = "capitata",
                                 "Brussels sprouts" = "gemmifera",
                                 "kohlrabi" = "gongylodes",
                                 "broccoli" = "italica")) %>% 
  left_join(., missingness) %>% 
  left_join(., sample_info) %>% 
  mutate(source = ifelse(id %in% genebank_ids,
                         paste0("Genebank"),
                         ifelse(grepl("Brassica_", population), 
                                paste0("wild"),
                                paste0(Type))),
         source = ifelse(source == "Germplasm",
                         "Genebank", 
                         paste0(source))) %>% 
  na.omit(cols = c("source")) %>% 
ggplot(., aes(LD1, LD2, color = population,
              shape = source)) + 
  geom_point(size = 3,
             alpha = 0.8) +
  scale_color_manual(values = cols) + 
  theme_minimal() # +
  facet_wrap(~population) 
  
  
  