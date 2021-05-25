# GWAS_Ole
GWAS of eye color

```r
library(dplyr)
library(ggplot2)
library(qqman)
library(SNPRelate)
library(gridExtra)
```
Making a binary phenotype file (brown vs other).

```r
nonbinary = read.table('eye_color.txt', header = F)
fam = read.table('eye_color.fam')
new.levels = c(2, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2)
FID = fam$V1
IID = nonbinary[,1]
phenotype = new.levels[factor(nonbinary$V2)]
write.table(data.frame(FID, IID, phenotype), file = 'binary_phenotype.txt', col.names = T, row.names = F, quote = F)
```

## QC 
#### Missing data

Testing for missing data

```bash 
plink --allow-no-sex --bfile eye_color --missing --out eye_color
```

This creates the files eye_color.imiss and eye_color.lmiss 

To analyse the heterozygosity of the sample --het is used

```bash
plink --allow-no-sex --bfile eye_color --het --out eye_color
```

The resulting .het file contains the observed number of homozygous genotypes as well as the number of non-missing genotypes per individual. With these values the observed heterozygosity rate per individual can be calculated using the this formula:

𝑜𝑏𝑠_{ℎ𝑒𝑡} = \frac{𝑁(𝑁𝑀) − 𝑂(𝐻𝑜𝑚)}{𝑁(𝑁𝑀)}

This is calculated in R and plotted with the proportion of missing SNPs per individuals (F_MISS in .imiss) as a function of the observed heterozygosity rate per individual.

```r
imiss = read.table('eye_color.imiss', header = T, sep="")
het = read.table('eye_color.het', header = T, sep="")

het$obs_het = (het$N.NM. - het$O.HOM.)/het$N.NM.
mean_het = mean(het$obs_het)
sd_het = sd(het$obs_het)

ggplot(data = het, aes(x = obs_het, y = imiss$F_MISS)) +
 geom_point(color='royalblue4', size=2, alpha=0.5) +
 labs(y = "Proportion of missing SNPs per individual", x = "Observed heterozygosity rate per individu
al") +
 geom_vline(xintercept = mean_het + 3*sd_het, color = "red") +
 geom_vline(xintercept = mean_het - 3*sd_het, color = "red") +
 theme_bw() +
 labs(title = "Heterozygosity rate versus proportion of missing SNPs") +
 theme(plot.title = element_text(hjust = 0.5))
```

The vertical lines indicate heterozygosity rates 3 standard deviations above and below the mean.
To filter out samples that are outside these lines, a file containing the FID and IID of these individuals is generated.

```r
right_tail = mean_het + 3*sd(het$obs_het)
left_tail = mean_het - 3*sd(het$obs_het)
filtering = cbind(imiss, het)
outlier_ind = subset(filtering, filtering$obs_het > right_tail | filtering$obs_het < left_tail)
write.table(outlier_ind[,c(1,2)], 'wrong_het_missing_values.txt', col.names = FALSE, row.names = FALSE)
```

This is then removed from the sample using plink

```bash
plink --allow-no-sex --bfile eye_color --remove wrong_het_missing_values.txt --pheno binary_phenotype.txt --make-bed --out eye_color_het
```
13 individuals were filtered out by doing this leaving us with 1274.

#### Relatedness

The remaining individuals in the sample are then analysed for relatedness, which is done by first
using the plink option –indep-pairwise. This creates a prune.in file, which is used to compute the IBD
matrix

```bash
plink --allow-no-sex --bfile eye_color_het --indep-pairwise 500kb 5 0.2 --out eye_color_het
plink --allow-no-sex --bfile eye_color_het --extract eye_color_het.prune.in --genome --min 0.185 --out eye_color_het
```
To remove the individuals with an IBD larger than 0.185, a .txt file is created in R.

```r
ibd = read.table('eye_color_het.genome', header = TRUE)
members = ibd[,1]
members = unique(members)
write.table(cbind(members,members), file = 'wrong_ibd.txt', col.names = F, row.names = F)
```
The individuals listed in wrong_ibd.txt are then removed using plink.

```bash
plink --allow-no-sex --bfile eye_color_het --remove wrong_ibd.txt --make-bed --out eye_color_het_ibd
```
1260 people are left in the data.

#### Variant filtering

Recomputing the missing data for each variant.

```bash
plink --allow-no-sex --bfile eye_color_het_ibd --missing --out eye_color_het_ibd
```

With the phenotype file, the Fisher’s exact test is performed on the case/control missing call counts

```bash
plink --allow-no-sex --bfile eye_color_het_ibd --pheno binary_phenotype.txt --test-missing --out eye_color_het_ibd
```

```r
test_missing = read.table('eye_color_het_ibd.missing', header = TRUE)
fail_diffmiss_qc = test_missing[test_missing$P < 10e-5, 2]
write.table(fail_diffmiss_qc, file = 'fail-diffmiss-qc.txt', row.names = F, col.names = F)
```

These variants are then filtered out along with any other that has a missing genotype rate larger than
0.05, deviates from Hardy-Weinberg equilibrium or has a minor allele frequency of less than 0.01.

```bash
plink --allow-no-sex --bfile eye_color_het_ibd --exclude fail-diffmiss-qc.txt --pheno binary_phenotype.txt --geno 0.5 --hwe 0.00001 --maf 0.01 --make-bed --out eye_color_het_ibd_var
```
This leaves 783902 variants.

### PCA
```r
snpgdsBED2GDS("eye_color_het_ibd_var.bed","eye_color_het_ibd_var.fam","eye_color_het_ibd_var.bim", "eye_color.gds")
```

Pruning the data for PCA

```bash
plink --allow-no-sex --bfile eye_color_het_ibd_var --indep-pairwise 500kb 5 0.2 --out eye_color_het_ibd_var
plink --allow-no-sex --bfile eye_color_het_ibd_var --pheno binary_phenotype.txt --extract eye_color_het_ibd_var.prune.in --pca 50 --out eye_color_het_ibd_var
```

```r
n_pcs <- min(1260, 783902)

genofile <- snpgdsOpen("eye_color.gds",  FALSE, TRUE, TRUE)
pca <- snpgdsPCA(genofile,eigen.cnt=n_pcs)
summary(pca)

eigenvectors <- as.data.frame(pca$eigenvect)
colnames(eigenvectors) <- as.vector(sprintf("PC%s", seq(1:nrow(pca$eigenvect))))
pca$sample.id <- sub("_eye_color_het_ibd_var_piece_dedup", "", pca$sample.id)
pca_percent <- pca$varprop*100

simpler.levels = c("Brown", "Other", "Other", "Other", "Other", "Other", "Brown", "Other", "Brown", "Other
", "Other", "Brown")

write.table(data.frame(FID, IID, eye_color = simpler.levels[factor(nonbinary$V2)]), file = 'simpler_eye_color.txt', col.names = T, row.names = F, quote = F)
info <- read.table("simpler_eye_color.txt", header = T)
eigenvectors$Phenotype <- info[match(pca$sample.id, info$IID),]$eye_color 

PC1_var <- pca$varprop[1]*100
PC2_var <- pca$varprop[2]*100

ggplot(data = eigenvectors, aes(x = PC1, y = PC2, col = Phenotype)) +
 geom_point(size=2,alpha=0.5) +
 scale_color_manual(breaks = c("Brown", "Other"), values=c("brown", "blue")) +
 xlab("PC1 (1.5%)") + ylab("PC2 (0.8%)") +
 labs(title = "PC1 x PC2") +
 theme(plot.title = element_text(hjust = 0.5))

ggplot(data = eigenvectors, aes(x = PC2, y = PC3, col = Phenotype)) +
    geom_point(size=2,alpha=0.5) +
    scale_color_manual(breaks = c("Brown", "Other"), values=c("brown", "blue")) +
    xlab("PC2 (0.8%)") + ylab("PC3 (0.3%)") +
    labs(title = "PC2 x PC3") +
    theme(plot.title = element_text(hjust = 0.5))
```

## Test for association

Fisher's exact test

```bash
plink --allow-no-sex --bfile eye_color_het_ibd_var --pheno binary_phenotype.txt --fisher --out eye_color_het_ibd_var_brown
```

p values plotted in Manhattan and QQ plot with Bonferroni corrected significance level

```r
association_test = read.table("eye_color_het_ibd_var_brown.assoc.fisher", header = T)
corrected_level <- 0.05/length(association_test$SNP)
manhattan(association_test, suggestiveline = FALSE, main = "Manhattan Plot", genomewideline = -log10(corrected_level), annotatePval = corrected_level)

r <- gcontrol2(association_test$P, col="black")
```

```r
(inflationfactor <- r$lambda)

association_test$GC  <- association_test$P/inflationfactor

gcontrol2(association_test$GC, col="black", main = "Corrected QQ plot")

manhattan(association_test, p = "GC", suggestiveline = FALSE, main = "Manhattan Plot with genomic contr
ol", genomewideline = -log10(corrected_level), annotatePval = corrected_level)
```

#### PC adjusting

First 20

```bash
plink --allow-no-sex --bfile eye_color_het_ibd_var --pheno binary_phenotype_brown.txt --logistic --covar eye_color_het_ibd_var.eigenvec --covar-number 1-20 --out eye_color_het_ibd_var_brown_20
```

QQplot and Manhattanplot

```r
log_regres_20 <- read.table("eye_color_het_ibd_var_brown_20.assoc.logistic", header = T)
only_snp_20 <- na.omit(subset(log_regres_20, TEST=="ADD"))
manhattan(only_snp_20, suggestiveline = FALSE, main = "Manhattan Plot", genomewideline = -log10(0.05/length(only_snp_20$SNP)), annotatePval = 0.05/length(only_snp_20$SNP))

cov_20 <- gcontrol2(only_snp_20$P, col="black", main = "QQ plot")

cov_20$lambda

log_regres_10 <- read.table("eye_color_het_ibd_var_brown_10.assoc.logistic", header = T)
only_snp_10 <- na.omit(subset(log_regres_10, TEST=="ADD"))
manhattan(only_snp_10, suggestiveline = FALSE, main = "Manhattan Plot", genomewideline = -log10(0.05/length(only_snp_10$SNP)), annotatePval = 0.05/length(only_snp_10$SNP))

cov_10 <- gcontrol2(only_snp_10$P, col="black", main = "QQ plot")

cov_10$lambda

log_regres_30 <- read.table("eye_color_het_ibd_var_brown_30.assoc.logistic", header = T)
only_snp_30 <- na.omit(subset(log_regres_30, TEST=="ADD"))
manhattan(only_snp_30, suggestiveline = FALSE, main = "Manhattan Plot", genomewideline = -log10(0.05/length(only_snp_30$SNP)), annotatePval = 0.05/length(only_snp_30$SNP))

cov_30 <- gcontrol2(only_snp_30$P, col="black", main = "QQ plot")

cov_30$lambda
```

Including 10 PCs seemed best as lambda was closest to 1 for the different PCs tested. 5, 15, 40 were also not better than 10.

Most significant SNP, 30 kb around it

```bash
plink --allow-no-sex --bfile eye_color_het_ibd_var --pheno binary_phenotype_brown.txt --recode A --snp rs1129038 --window 30 --out eye_color_het_ibd_var_brown
```

Genotype distributions

```r
d <- read.table("eye_color_het_ibd_var_brown.raw", header = T)
d <- d %>% mutate(pheno = ifelse(PHENOTYPE == 2, "Brown", "Other"))

snp0 <- d %>% filter(d$rs1129038_C == 0)
snp1 <- d %>% filter(d$rs1129038_C == 1)
snp2 <- d %>% filter(d$rs1129038_C == 2)

l <- ggplot(snp0, aes(x = pheno, fill = pheno)) +
 geom_bar(fill = c('sienna4', 'darkseagreen3')) +
 theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 11)) + guides(fill=FALSE) + ggtitle("Homozygous")

m <- ggplot(snp1, aes(x = pheno, fill = pheno)) +
 geom_bar(fill = c('sienna4', 'darkseagreen3')) +
 theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 11)) + guides(fill=FALSE) + ggtitle("Heterozygous")

n <- ggplot(snp2, aes(x = pheno, fill = pheno)) +
 geom_bar(fill = c('sienna4', 'darkseagreen3')) +
 theme(axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, size = 11)) + guides(fill=FALSE) + ggtitle("Homozygous")

grid.arrange(l, m, n, ncol=3, top = "Distribution of each genotype")
```

# Epistasis test

Doing an exhaustive epistasis test with the faster --fast-epistasis command.

```bash
plink --allow-no-sex --bfile eye_color_het_ibd_var --fast-epistasis --out eye_color_brown
```
Calculating the distance between each test pair and plotting the result.

```r
f <- read.table("eye_color_brown.epi.cc", header = T)

snp1 <- log_regres_10[match(f$SNP1, log_regres_10$SNP),]$BP
snp2 <- log_regres_10[match(f$SNP2, log_regres_10$SNP),]$BP
f$BP1 <- snp1
f$BP2 <- snp2
bonf <- 0.001/length(f$P)
f$Distance <- abs(snp2-snp1)

ggplot(f) +
    geom_point(aes(x = Distance, y = P, col=P < bonf)) +
    scale_color_manual(breaks = c("P > α", "P < α"), values = c('black', 'darkseagreen3'))

```

Finding the 100 SNP pairs with the smallest p-values.

```r
df <- f %>% 
    arrange(P) %>%
    top_n(-100) %>% 
    select(SNP1, SNP2)
```

Making a file for another epistasis test.
```r
write.table(df, "epi_x.set", append = FALSE, sep = " ", dec = ".", row.names = FALSE, col.names = FALSE)
```
Running the more precise --epistasis test on these 100 pairs

```bash
plink --allow-no-sex --bfile eye_color_het_ibd_var --epistasis --set-test --set epi.set --out eye_color_top100
```
```r
top100 <- read.table("eye_color_top100.epi.cc", header = T)

snp1 <- log_regres_10[match(top100$SNP1, log_regres_10$SNP),]$BP
snp2 <- log_regres_10[match(top100$SNP2, log_regres_10$SNP),]$BP
top100$BP1 <- snp1
top100$BP2 <- snp2
bonf <- 0.001/length(top100$P)
top100$Distance <- abs(snp2-snp1)

ggplot(top100) +
    geom_point(aes(x = Distance, y = P, col=P < bonf)) +
    scale_color_manual(breaks = c("P > α", "P < α"), values = c('black', 'darkseagreen3'))
   
```
```r
ggplot(top100, aes(x = SNP1, y = SNP2, fill = P)) +
    geom_tile()
```

Grouping the chromosomes together, calculating the mean p-value and plotting the result.
```r
top100_chr <- top100 %>%
    group_by(CHR1, CHR2) %>% 
    summarise(chr_mean = mean(P))

ggplot(top100_chr, aes(x = CHR1, y = CHR2, fill = chr_mean)) +
    geom_tile() +
    labs(fill='Mean p-value')
```
