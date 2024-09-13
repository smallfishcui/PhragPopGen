cd /cuiwang/nodelete/rawRAD/Dstat
ln -s /cuiwang/nodelete/rawRAD/populations/firstsnpperloci/populations.snps.vcf .
vcftools --vcf populations.snps.vcf --out ref_structure --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --remove-indv Y56 --recode
module load admixtools
sh convertVCFtoEigenstrat.sh --renameScaff
vcftools --gzvcf ref_structure.recode.renamedScaff.vcf.gz --out ref_structure --plink
convertf -p par.PED.EIGENSTRAT.ref_structure.recode
awk 'BEGIN{i=0}{i=i+1; print $1"\t"$2"\t"$3"\t"i"\t"$5"\t"$6}' ref_structure.snp > ref_structure.snp.tmp
mv ref_structure.snp.tmp ref_structure.snp
awk '{ print $1"\t"$2 }' ref_structure.ind > ref_structure.ind_bak
#add population label to the 3rd column
.libPaths(c("/projappl/project_2001259/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
#install.packages("remotes")
#devtools::install_github("uqrmaie1/admixtools")
library(admixtools)
library(tidyverse)
prefix ='./ref_structure'
my_f2_dir = './resadmixtools'
extract_f2(prefix, my_f2_dir)
f2_blocks = f2_from_precomp(my_f2_dir)
F3=f3(f2_blocks)
write.table(F3, file = "f3_for_all.txt", sep = "\t", row.names=FALSE, col.names=TRUE, quote=FALSE )
#f4
popA1='AU'
popB1=c('CN', 'Med','Usland','EU','Usnat')
popC1=c('CN', 'Med','Usland','EU','Usnat')
popD1=c('CN', 'Med','Usland','EU','Usnat')
#pops = c('AU','CN', 'Med','Usland','EU','Usnat')
qpstatAU=qpdstat(f2_blocks, popA1, popB1, popC1, popD1)
qpstat=qpdstat(f2_blocks, pops,pops,pops,pops)
write.table(qpstat, file = "qpstat.txt", sep = " ", row.names=FALSE, col.names=TRUE, quote=FALSE )

#plot F3 for one population;https://genoplot.com/discussions/topic/28207/r-scripts-for-admixtools-2/1
library(pheatmap)
library(colorspace) # for hex
library(vegan) # for reorder.hclust

p1="Usland"
p2=c('AU','CN', 'Med','EU','Usnat')
f3=f3(f2_blocks,p1,p2,p2)

t=f3%>%select(pop2,pop3,est)%>%pivot_wider(names_from=pop2,values_from=est)%>%data.frame(row.names=1)

diag(t)=NA # don’t show the f3 values of populations with themselves in order to not stretch the color scale

pheatmap(t,filename="Usland_f3.png", display_numbers=apply(t,1,function(x)sub("^0","",sprintf("%.3f",x))), legend=F,border_color=NA,cellwidth=16,cellheight=16,fontsize=8,fontsize_number=7,number_color=“black”, breaks=seq(.160,.190,.03/256), \ 
# use fixed limits for color scale \
colorRampPalette(hex(HSV(c(210,170,120,60,40,20,0),.5,1)))(256))
pheatmap(t,filename="Usland_f3.png",display_numbers=apply(t,1,function(x)sub("^0","",sprintf("%.3f",x))),legend=F,border_color=NA,cellwidth=16,cellheight=16,fontsize=8,fontsize_number=7)
#fst
fst(my_f2_dir)
