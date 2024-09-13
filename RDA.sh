#create genotype files
cd /scratch/project_2001259/cuiwang/nodelete/rawRAD/populations/firstsnpperloci
mkdir RDA
vcftools --vcf populations.snps.vcf --out RDA_allpure --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --remove-indv Y56 --remove-indv K13 --remove-indv Y43 --remove-indv Y48 --remove-indv Y49 --remove-indv Y55 --remove-indv Y44 --remove-indv Y36 --remove-indv F2 --remove-indv H2 --remove-indv Y3 --remove-indv K12 --plink
module load plink/1.90
#remove the three individuals that had wrong coordinates
plink --file RDA_allpure --recode A --noweb --remove mylist.txt #note here don’t use recode AD

awk 'NR==FNR{a[$1]=$0; next} {print a[$1]}' "RDA_sampleinfo.txt" "plink.nosex" | awk 'NF' > RDA_sampleinfo.sorted.txt

#puhti
module load r-env
start-r 
library(vegan)
library(psych)
library(adegenet)
GenData = read.PLINK("plink.raw", row.names=1)
ClimData=read.table(file="./RDA_sampleinfo.sorted.txt", header = T, dec = ",")
gen.imp <- apply(GenData, 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))

dim(gen.imp)
sum(is.na(gen.imp))
str(ClimData)
ClimData$Sample<-as.character(ClimData$Sample)
ClimData$Group<-as.factor(ClimData$Group)
ClimData$Ploidy<-as.factor(ClimData$Ploidy)
ClimData$lat<-as.numeric(ClimData$Latitude)
ClimData$long<-as.numeric(ClimData$Longititude)
#the order of the files are important
identical(rownames(gen.imp), ClimData[,1])

#RDA1 <- rda(gen.imp ~ lat + long + Condition(Ploidy), ClimData), constrait 9.34, conditional 15.09
#or RDA1 <- rda(gen.imp ~ lat + long + Condition(Group), ClimData), constraint 3.98, conditional 44.97
LL=cbind(ClimData$lat,ClimData$long)
#check the group proportion, result is constraint 34.11, conditional 14.84
RDA1 <- rda(gen.imp ~ ClimData$Group + Condition(LL))
# check the ploidy proportion,result is Constrained 9.59, Conditional 14.84
RDA1 <- rda(gen.imp ~ ClimData$Ploidy + Condition(LL))
#check the ploidy without group, result is Conditional 44.97; Constrained 2.445.
RDA1 <- rda(gen.imp ~ ClimData$Ploidy + Condition(ClimData$Group))
#check the ploidy without group, result is Conditional 15.09; Constrained 32.32.
RDA1 <- rda(gen.imp ~ ClimData$Group + Condition(ClimData$Ploidy))

#use eiher group or ploidy as a covariate
RsquareAdj(RDA1)
summary(eigenvals(RDA1, model = "constrained"))
summary(eigenvals(RDA1, model = "unconstrained"))
screeplot(RDA1)
signif.axis<-anova.cca(RDA1, permutations = how(nperm=100), by="terms", parallel=8)
signif.axis
vif.cca(RDA1)
#levels(ClimData$Ploidy)<-c("4x","6x","8x")
levels(ClimData$Group)<-c("AU", "CN", "EU", "Med", "Usland", "Usnat")
#ploidy<-ClimData$Ploidy
group<-ClimData$Group
#bg <- c("dark red","dark blue","orange")
bg <- c("#F1605DFF","#FB9A06FF","#35B779FF","#00336FFF","#7301A8FF","#CD4071FF") 
# axes 1 & 2
pdf("RDA1_2.pdf")
plot(RDA1, type="n", scaling=3)
#points(RDA1, display="species", pch=20, cex=0.7, col="gray32", scaling=3)# the SNPs
points(RDA1, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=bg[group]) # the samples
text(RDA1, scaling=3, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottomright", legend=levels(group), bty="n", col="gray32", pch=21, cex=1, pt.bg=bg)
#the proportion
R.sum <- summary(RDA1)
R.sum$cont   # Prints the "Importance of components" table
R.sum$cont$importance[2, "RDA1"]#7.2
R.sum$cont$importance[2, "RDA2"]#3.7
#install.packages('ggord'), https://fukamilab.github.io/BIO202/06-B-constrained-ordination.html
library(ggord)
library(ggplot2)
pdf("shaded.pdf")
ggord(RDA1, ClimData$Group) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
dev.off()

