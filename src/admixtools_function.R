################################################################################
## D statistics using ADMIXTOOLS & R/admixr
## Author: Sarah Turner-Hissong
## Date: 26 November 2019
## follows tutorial at https://bodkan.net/admixr/articles/tutorial.html
################################################################################

# This script computes D statistics using the ADMIXTOOLS software
# (https://github.com/DReichLab/AdmixTools) which is called by R/admixr

library(admixr)
library(tidyverse)

################################################################################
## load data (EIGENSTRAT format)
## ind, snp, and geno files
################################################################################

## command line code to create input file
## plink --vcf ../../data/processed/filtered_snps/oleracea_filtered.vcf.gz \
## --recode \
## --out oleracea_filtered \
## --const-fid \
## --allow-extra-chr
## sed "s/C//g" data/oleracea_filtered.map > data/oleracea.map
## ~/software/EIG-7.2.1/src/convertf -p convertf_parfile

# read in SNP data
snps <- eigenstrat(prefix = "data/oleracea")

# define population 3
pop3 <- "palmifolia"

# define populations
pops <- c("alboglabra", "botrytis", "capitata", "sabauda", "costata",
          "gemmifera", "gongylodes", "italica", "oleracea", "palmifolia",
          "medullosa", "ramosa", "sabellica", "viridis")

pops2 <- pops[!pops %in% pop3]

pops <- c("alboglabra", "botrytis", "sabauda", "costata",
          "gemmifera", "gongylodes", "italica", "oleracea",
          "medullosa", "ramosa", "viridis")


pops2 <- c("alboglabra", "botrytis", "capitata", "sabauda", "costata",
          "gemmifera", "gongylodes", "italica", "oleracea",
          "medullosa", "ramosa", "viridis")

# calculate D statistics (ABBA-BABA)

D_result1 <- d(W = pops,
            X = "capitata",
            Y = "palmifolia",
            Z = "Brassica_rupestris",
            data = snps)

D_result2 <- d(W = pops2,
            X = "sabellica",
            Y = "palmifolia",
            Z = "Brassica_rupestris",
            data = snps)

D_result3 <- d(W = pops2,
            X = "palmifolia",
            Y = "sabellica",
            Z = "Brassica_rupestris",
            data = snps)

# D_result2 <- d(W = pops,
#             X = "gongylodes",
#             Y = "sabellica",
#             Z = "Brassica_rupestris",
#             data = snps)

write.csv(D_result1, "results/capitata_palmifolia.csv")
write.csv(D_result2, "results/sabellica_palmifolia.csv")
write.csv(D_result3, "results/palmifolia_sabellica.csv")


# D_result

Dplot1 <- ggplot(D_result1, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Population")

ggsave("results/capitata_palmifolia.png", height = 4, width = 6)

Dplot2 <- ggplot(D_result2, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Population")

ggsave("results/sabellica_palmifolia.png", height = 4, width = 6)

Dplot2 <- ggplot(D_result3, aes(fct_reorder(W, D), D, color = abs(Zscore) > 2)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_errorbar(aes(ymin = D - 2 * stderr, ymax = D + 2 * stderr)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Population")

ggsave("results/palmifolia_sabellica.png", height = 4, width = 6)

# f4 statistics
# departure from 0 indicates evidence of gene flow

# f4_result <- f4(W = pops2,
#               X = "italica",
#               Y = pop3,
#               Z = "Brassica_rupestris",
#               data = snps)
#
# f4_result
#
# f4plot <- ggplot(f4_result, aes(fct_reorder(W, f4), f4, color = abs(Zscore) > 2)) +
#   geom_point() +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_errorbar(aes(ymin = f4 - 2 * stderr, ymax = f4 + 2 * stderr)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
#   xlab("Population")
#
# ggsave("results/f4plot.png", height = 4, width = 6)
#
# # f4 ratio statistic
# # proportion of shared ancestry
#
# f4_ratio <- f4ratio(X = pops2,
#               A = "italica",
#               B = "palmifolia",
#               C = "Brassica_insularis",
#               O = "Brassica_rupestris",
#               data = snps)
#
# f4ratioplot <- ggplot(f4_ratio, aes(fct_reorder(X, alpha), alpha, color = abs(Zscore) > 2)) +
#   geom_point() +
#   geom_hline(yintercept = 0, linetype = 2) +
#   geom_errorbar(aes(ymin = alpha - 2 * stderr, ymax = alpha + 2 * stderr)) +
#   xlab("Population")
#
# ggsave("results/f4_ratio_plot.png", height = 4, width = 6)
#
# # f3 statistic
#
# f3result <- f3(A = pops,
#             B = pops,
#             C = "Brassica_rupestris",
#             data = snps)
#
# f3result
#
# f3plot <- f3result %>%
#   filter(A != B) %>%
#   # mutate(A = factor(A, levels = ordered),
#   #        B = factor(B, levels = ordered)) %>%
#   ggplot(aes(A, B)) + geom_tile(aes(fill = f3))
#
# ggsave("results/f3plot.png", height = 5, width = 5)
