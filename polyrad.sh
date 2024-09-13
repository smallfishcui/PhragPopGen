cd /cuiwang/nodelete3/rawRAD/populations/firstsnpperloci
bcftools sort -m 76000M -o populations.snps.sorted.vcf.gz -Oz -T ./ populations.snps.vcf
bcftools view -S sample_noOG_USland.txt -Oz -T ./ -o populations.snps.sorted.noOGUSland.vcf.gz populations.snps.sorted.vcf.gz
cd /cuiwang/nodelete3/RADoctoploid

module load r-env
start-r
.libPaths(c("/projappl/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
library(polyRAD)
library(adegenet)
library(polyRADtutorials)
library(ggplot2)
library(maps)
library(pegas)
library(PBSmapping)
library(spdep)
library(polysat)
library(tidyr)
library(umap)

#info<- read.delim("/cuiwang/nodelete3/RDA.txt", sep="\t", header=T)
info<- read.delim("./info_noOGUSland.txt", sep="\t", header=T)
# Make a vector of taxa ploidies from our metadata spreadsheet
tp <- info$Ploidy
names(tp) <- info$Sample
mydata<-VCF2RADdata("/cuiwang/nodelete3/rawRAD/populations/firstsnpperloci/populations.snps.sorted.noOGUSland.vcf.gz",phaseSNPs=FALSE, refgenome ="/cuiwang/nodelete3/rawRAD/PaEUcurated.fasta", min.ind.with.reads=44, min.ind.with.minor.allele=1,taxaPloidy = tp, expectedLoci=1524048)
#mydata<-VCF2RADdata("/cuiwang/nodelete3/rawRAD/populations/firstsnpperloci/populations.snps.sorted.vcf.gz",phaseSNPs=FALSE, refgenome ="/cuiwang/nodelete3/rawRAD/PaEUcurated.fasta", min.ind.with.reads=44, min.ind.with.minor.allele=1,taxaPloidy = tp, expectedLoci=1524048)

mydata <- AddPCA(mydata)
hh <- HindHe(mydata)
info$HindHe <- rowMeans(hh, na.rm = TRUE)[info$Sample]
info$Depth <- rowSums(mydata$locDepth)[info$Sample]
info$Fis <- NA
for (i in 1:nrow(info)) {
  info$Fis[i] <- InbreedingFromHindHe(info$HindHe[i], info$Ploidy[i])
}

pdf("depth_he.pdf")
ggplot(info, aes(x = Depth, y = HindHe)) +
  geom_point() +
  facet_wrap(~ Ploidy) +
  ggtitle("Read depth and Hind/He across individuals") +
  scale_x_log10()
dev.off()

mydata_tetra <- SubsetByPloidy(mydata, ploidies = mydata$possiblePloidies[4])
od <- TestOverdispersion(mydata_tetra, to_test = 8:20) #optimal=9

alfreq <- colMeans(mydata$depthRatio, na.rm = TRUE)
theseloci2x <- GetLoci(mydata)[mydata$alleles2loc[alfreq >= 0.05 & alfreq < 0.5]]
theseloci2x <- unique(theseloci2x)
hh2x_05 <- colMeans(hh[, theseloci2x], na.rm = TRUE)

info$HindHe <- rowMeans(hh[, theseloci2x], na.rm = TRUE)[info$Sample]
info$Fis <- NA
for (i in 1:nrow(info)) {
  info$Fis[i] <- InbreedingFromHindHe(info$HindHe[i], info$Ploidy[i])
}

alfreq <- colMeans(mydata$depthRatio, na.rm = TRUE)
theseloci <- GetLoci(mydata)[mydata$alleles2loc[alfreq >= 0.05  & alfreq < 0.5]]
theseloci <- unique(theseloci)
hh4x_05 <- colMeans(hh[mydata$taxaPloidy == 4, theseloci], na.rm = TRUE)
pdf("tetraploid_hh.pdf")
hist(hh4x_05, breaks = 20, xlab = "Hind/He", main = "Hind/He in tetraploids, MAF >= 0.05")
dev.off()

mydata_tetraPopStruct <- IteratePopStruct(mydata_tetra, overdispersion = od$optimal)

hhByLoc <- colMeans(hh, na.rm = TRUE)
set.seed(528)
ExpectedHindHe(mydata_tetra, inbreeding = 0.65, overdispersion = od$optimal)
thresh1 <- 0
thresh2 <- 1.06
keeploci <- names(hhByLoc)[hhByLoc > thresh1 & hhByLoc < thresh2]

mydata2 <- SubsetByLocus(mydata, keeploci)
mydata2

set.seed(326)
mydata2PopStruct <- IteratePopStruct(mydata2, overdispersion = od$optimal)
identical(rownames(mydata2PopStruct$PCA), info$Sample)
info<- cbind(info, mydata2PopStruct$PCA)
matrixPopStruct <- GetWeightedMeanGenotypes(mydata2PopStruct)
clust1 <- find.clusters(matrixPopStruct, n.pca = nrow(matrixPopStruct))
dapc_PopStruct <- dapc(matrixPopStruct, clust1$grp,
                       n.pca = 200, n.da = 6)
info$DAPC_PopStruct <- dapc_PopStruct$assign
table(info$DAPC_PopStruct, info$Ploidy)
pdf("hhDAPC_PopStruct.pdf")
ggplot(info, aes(x = DAPC_PopStruct, y = HindHe, fill = DAPC_PopStruct)) +
  geom_boxplot()
dev.off()
pdf("UMAP.pdf")
ggplot(info, aes(x = PC1, y = PC2, color = Group, shape = as.factor(Ploidy))) +
  geom_point()
dev.off()

# Subset the data based on the given criteria
subset_data <- subset(info, 
                      (Group == "AU" & Ploidy == 4) |
                      (Group == "Usland" & Ploidy == 3) |
                      (Group != "AU" & Group != "Usland" & Ploidy == 2))

stats <- aggregate(Fis ~ Group, data = subset_data, 
                     FUN = function(x) c(mean = mean(x), sd = sd(x)))
  