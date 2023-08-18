################################################################################
## Diagnostic plots for variants
## Taken from the excellent tutorial at https://evodify.com/gatk-in-non-model-organism/
## (with some minor modifications)
################################################################################

library(tidyverse)
library(cowplot)

VCFsnps <- read.csv('reports/filtering_bpres/all_samps.raw.snps.table',
                    header = T,
                    na.strings=c("","NA"),
                    sep = "\t")

dim(VCFsnps)

VCFinvariant <- read.csv('reports/filtering_bpres/all_samps.raw.invariant.table',
                         header = T,
                         na.strings=c("","NA"),
                         sep = "\t")

dim(VCFinvariant)

VCF <- bind_rows(list(variant = VCFsnps,
                      invariant = VCFinvariant),
                 .id = "class")

dim(VCF)

xmax <- mean(VCF$DP)*3

# colors
snps <- '#A9E2E4'
invariant <- '#F4CCCA'

DP <- ggplot(VCF, aes(x=DP, fill=class)) +
  geom_density() +
  geom_vline(xintercept=c(10,xmax)) +
  xlim(c(0, 6000)) +
  theme_bw()

QD <- ggplot(VCF, aes(x=QD, fill=class)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept=2, size=0.7) +
  theme_bw()

FS <- ggplot(VCF, aes(x=FS, fill=class)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept=c(60, 200), size=0.7) +
  ylim(0,0.1) +
  theme_bw()

MQ <- ggplot(VCF, aes(x=MQ, fill=class)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept=40, size=0.7) +
  theme_bw()

MQRankSum <- ggplot(VCF, aes(x=MQRankSum, fill=class)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept=-20, size=0.7) +
  theme_bw()

SOR <- ggplot(VCF, aes(x=SOR, fill=class)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept=c(4, 10), size=1, colour = c(snps, invariant)) +
  theme_bw()

ReadPosRankSum <- ggplot(VCF, aes(x=ReadPosRankSum, fill=class)) +
  geom_density(alpha=.3) +
  geom_vline(xintercept=c(-10,10), size=1) +
  xlim(-30, 30) +
  theme_bw()

plot_grid(QD, DP, FS, MQ, MQRankSum, SOR, ReadPosRankSum, nrow=4)

ggsave("reports/filtering_bpres/unfiltered_diagnostics.png")
