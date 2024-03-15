library(adegenet)
library(pegas)
library(vcfR)
library(tidyverse)
library(readxl)

# samples <- read.table("~/Documents/bo-demography/models/pca/pca_results.txt", header = TRUE)
samples <- read.table("~/Documents/bo-demography/models/pca/tassel_mds.txt", header = TRUE, skip = 2, stringsAsFactors = FALSE)

str(samples)

samples <- samples %>%
  separate(Taxa, into = c("id", "cell", "run"), sep = "_", fill = "right", remove = FALSE) %>%
  mutate(id = ifelse(!grepl("Bo", id), paste0("SamC_", cell), id)) %>%
  select(Taxa, id)

sample_ids <- read_xlsx("~/Documents/brassica-oleracea/Bo_sample_ids.xlsx") %>%
  mutate(id = paste0("Bo", id, rep)) %>%
  select(id, "Genus", "Species", "Variety/Form") %>%
  mutate(morph = ifelse(is.na(`Variety/Form`), paste0("Brassica_", Species), `Variety/Form`)) %>%
  select(id, morph) %>%
  mutate(morph = fct_recode(morph, viridis = "virdis"))

sample_ids2 <- data.frame(id = paste0("SamC_", str_pad(1:119, width = 3, pad = "0")),
                          morph = c(rep("capitata", 45), 
                                    rep("gongylodes", 19),
                                    rep("botrytis", 20),
                                    rep("italica", 23),
                                    rep("alboglabra", 4),
                                    rep("gemmifera", 2),
                                    rep("acephala", 2),
                                    rep("sabellica", 2),
                                    rep("oleracea", 2)))

sample_ids <- rbind(sample_ids, sample_ids2) %>%
  left_join(samples, .) %>%
  mutate(morph = fct_recode(morph, palmifolia = "longata"))

sample_ids <- sample_ids %>%
  mutate(group = ifelse(morph %in% c("capitata", "costata", "medullosa", "sabauda"), "cabbage", 
                        ifelse(morph %in% c("longata", "ramosa", "sabellica", "viridis", "acephala"), "kale",
                               ifelse(morph %in% c("alboglabra", "botrytis", "italica"), "floral",
                                      ifelse(morph %in% "gongylodes", "kohlrabi", 
                                             ifelse(morph %in% "gemmifera", "Brussels sprouts", paste(morph)))))),
         morph = ifelse(Taxa %in% "Bo167A_S70_L001", "botrytis", paste(morph)))

sample_ids %>%
  select(morph, Taxa) %>%
  write.table("~/Documents/bo-demography/models/smc/population_ids.txt",
              row.names = FALSE,
              quote = FALSE)

oleracea <- read.vcfR("~/Documents/bo-demography/models/pca/oleracea_thinned.recode.vcf")

oleracea.genlight <- vcfR2genlight(oleracea)

# plot SFS
mySum <- glSum(oleracea.genlight, alleleAsUnit = TRUE)
barplot(table(mySum), col="blue", space=0, xlab="Allele counts")

pca.1 <- glPca(oleracea.genlight, nf=300) 


# pca
pca1 <- readRDS("~/Documents/bo-demography/models/pca/pca.rds")

str(pca1)
pop(oleracea.genlight) <- sample_ids$morph 

g1 <- s.class(pca1$scores, pop(oleracea.genlight), xax=1, yax=2)

pcs <- as.data.frame(pca1$scores) %>%
  cbind(., sample_ids) %>%
  mutate(type = ifelse(grepl("Brassica|oleracea", morph), "wild", "cultivated"),
         wild = ifelse(grepl("Brassica", morph), paste(morph), "cultivated"),
         morph = fct_recode(morph, "Chinese kale" = "alboglabra",
                            "cauliflower" = "botrytis",
                            "cabbage" = "capitata",
                            "tronchuda cabbage" = "costata",
                            "Brussels sprouts" = "gemmifera",
                            "kohlrabi" = "gongylodes",
                            "broccoli" = "italica",
                            "marrow cabbage" = "medullosa",
                            "Tuscan kale" = "palmifolia",
                            "perpetual kale" = "ramosa",
                            "savoy cabbage" = "sabauda",
                            "curly kale" = "sabellica",
                            "collards" = "viridis",
                            "cabbage" = "acephala"),
         cultivated = ifelse(grepl("Brassica|oleracea", morph), "wild", paste(morph)))

ggplot(pcs, aes(PC2, PC3, shape = type)) +
  geom_point(data = subset(pcs, type %in% "wild"),
             aes(PC2, PC3), colour = "gray",
             size = 3, alpha = 0.8, shape = "triangle") + 
  geom_point(data = subset(pcs, type %in% "cultivated"),
             aes(PC2, PC3, colour = morph),
             size = 3, alpha = 0.8) + 
  # scale_colour_brewer(palette = "Set1") + 
  theme_minimal() +
  theme(legend.title = element_blank()) +
  # guides(colour = guide_legend(override.aes = list(size = 1.5))) + 
  # xlab("PC1 (10.6% variance explained)") + 
  xlab("PC2 (4.6% variance explained)") + 
  ylab("PC3 (3% variance explained)")

ggsave("~/Desktop/pc2_pc3.png", height = 4.5, width = 6.5)

pca1$eig[1]/sum(pca1$eig) # proportion of variation explained by 1st axis 
pca1$eig[2]/sum(pca1$eig)
pca1$eig[3]/sum(pca1$eig)



oleracea.genind <- vcfR2genind(oleracea)

X <- scaleGen(oleracea.genind, NA.method = "mean")

dim(X)
X[1:5,1:5]

pca2 <- prcomp(X)

library(ggfortify)
autoplot(pca2, x = 2, y = 3)

grp <- find.clusters(oleracea.genlight, max.n.clust = 40)

pca1 <- dapc(oleracea.genlight, grp$grp)

table.value(table(sample_ids$morph, grp$grp), col.lab=paste("inf", 1:23),
            row.lab=paste("ori", 1:23))


pca2 <- glPca(oleracea.genlight, nf = 50)

scatter(pca2, 1:2)

scatter(pca2, 2:3)

library(ape)

mycol <- colorplot(pca2$scores[,1:2], pca2$scores[,1:2], transp = TRUE, cex = 2)
colorplot(pca2$scores[,2:3], pca2$scores[,2:3], transp = TRUE, cex = 2)

s.class(pca2$loadings, fac = as.factor(sample_ids$morph), col = funky(15))

tree <- nj(dist(as.matrix(oleracea.genlight)))

library(ggtree)

ggtree(tree, layout = "unrooted") %<+%
  sample_ids +
  geom_tiplab(aes(color = morph), size = 2) +
  theme(legend.position = "right") +
  geom_label2(aes(subset = !isTip, label = node), size = 1)


drop <- sample_ids[grepl("Brassica", sample_ids$morph), "Taxa"]
drop <- drop[-27]


pruned.tree <- drop.tip(tree, tree$tip.label[match(drop, tree$tip.label)])

tree2 <- root(pruned.tree, out = 146)

ggtree(tree2, branch.length = "none", layout = "circular") %<+%
  sample_ids +
  geom_tiplab2(aes(color = morph), size = 2) +
  # theme(legend.position = "right") +
  geom_hilight(node = 282, fill = "darkgreen", alpha = 0.6) + 
  geom_hilight(node = 361, fill = "orange", alpha = 0.6) + 
  geom_hilight(node = 396, fill = "steelblue", alpha = 0.6) +
  geom_hilight(node = 273, fill = "lightgreen", alpha = 0.6) + 
  geom_hilight(node = 437, fill = "darkorange", alpha = 0.6) +
  geom_hilight(node = 427, fill = "yellow", alpha = 0.6) +
  geom_hilight(node = 450, fill = "blue", alpha = 0.6) +
  geom_hilight(node = 349, fill = "purple", alpha = 0.6) +
  geom_hilight(node = 512, fill = "firebrick", alpha = 0.6) +
  geom_hilight(node = 504, fill = "red", alpha = 0.6) +
  geom_hilight(node = 489, fill = "pink", alpha = 0.6) 
  
#   geom_label2(aes(subset = !isTip, label = node), size = 3)


ggsave("~/Desktop/rooted_circular.pdf", height = 15, width = 15)

plot(tree, typ="unrooted", 
     show.tip = FALSE, 
     edge.width = 2, 
     edge.color = sample_ids$morph, 
     col = mycol)

tiplabels(pch = 20, col = sample_ids$morph)
tiplabels(pch = 20, col = mycol, cex = 2)

plot(tree, type = "unrooted")

# inbreeding
temp <- inbreeding(oleracea.genlight, N = 2)
Fbar <- sapply(temp, mean)
Fbar

Fbar %>%
  as.data.frame(Fbar) %>%
  ggplot(., aes(Fbar)) +
  geom_histogram()

# pca
sum(is.na(genlight$tab))

X <- scaleGen(genlight, NA.method = "mean")

dim(X)
X[1:5,1:5]

pca1 <- dudi.pca(X, cent = FALSE, scale = FALSE, scannf = FALSE, nf = 50)

barplot(pca1$eig[1:50])

pca1

s.label(pca1$li)

ggplot(pca1$li, aes(Axis1, Axis2)) +
  geom_point()



