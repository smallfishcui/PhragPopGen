module load biokit
bcftools view -s C2,C3,C6,C7,H1,Y1,Y25,Y39,Y51_1,Y53,Y6,Y9 populations.snps.vcf > AU.vcf
bcftools view -s Y32,Y38,Y40,Y41,Y45,Y46,Y47,Y57 populations.snps.vcf > CN.vcf
bcftools view -s H6,K2,A10,A2,A3,A5,A7,A8,B2,C10,C11,C12,C13,C9,D11,D12,D13,D4,D6,D7,D9,E3,F13,Y10,Y11,Y13,Y14,Y15,Y16,Y19,Y20,Y22,Y23,Y24,Y31,Y4 populations.snps.vcf > EU.vcf
bcftools view -s I12,A9,D5,E9,E13,F9,H5,Y18,Y21,Y26,Y28,Y34 populations.snps.vcf > Med.vcf
bcftools view -s Y7,Y12,Y30,Y33,Y35,Y37 populations.snps.vcf > Usland.vcf
bcftools view -s Y17_1,E12,Y27 populations.snps.vcf > Usnat.vcf

#
#vk tajima 10,000 10,000 Usnat.vcf > Usnat.tajima
#use ucftools
vcftools --vcf CN.vcf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --TajimaD 10000 --out CN 
vcftools --vcf Med.vcf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --TajimaD 10000 --out Med
vcftools --vcf Usnat.vcf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --TajimaD 10000 --out Usnat
vcftools --vcf Usland.vcf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --TajimaD 10000 --out Usland
vcftools --vcf EU.vcf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --TajimaD 10000 --out EU
vcftools --vcf AU.vcf --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --TajimaD 10000 --out AU
paste AU.Tajima.D CN1.Tajima.D EU1.Tajima.D Med1.Tajima.D Usland1.Tajima.D Usnat1.Tajima.D > all.Tajima.D

#modify the txt, use ggplot2 to plot the violin statistic
library(ggplot2)
a<-read.table("plotTajima.txt",header=TRUE)
a$Lineage<-as.factor(a$Lineage)
pdf("Tajima_lineage.pdf")
p<-ggplot(a, aes(x=Lineage, y=TajimaD, fill=Lineage)) +
  geom_violin(trim=FALSE)
p + stat_summary(fun.data=mean_sdl, mult=1, 
                 geom="pointrange", color="red")
dev.off()

library(psych)
au<-read.table("plotAU.txt",header=TRUE)
describe(au)
EU<-read.table("plotEU.txt",header=TRUE)
describe(EU)
CN<-read.table("plotCN.txt",header=TRUE)
describe(CN)
Med<-read.table("plotMed.txt",header=TRUE)
describe(Med)
Usland<-read.table("plotUsland.txt",header=TRUE)
describe(Usland)
Usnat<-read.table("plotUsnat.txt",header=TRUE)
describe(Usnat)