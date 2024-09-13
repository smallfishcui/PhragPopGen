cd /cuiwang/nodelete/rawRAD/process_radtags
#discard reads that have multiple alignments
for i in `ls *.bam | sed 's/\.bam//'`; do samtools view -@40 -F 2048 -bo $i".filtered.bam" $i".bam"; done
#move the unfiltered bam file to another folder
awk '{ print $1"\t"$2 }' /cuiwang/nodelete/popmap_denovo20190605_sorted.txt > ./popmap20240123.txt
ref_map.pl --samples /cuiwang/nodelete/rawRAD/process_radtags/primaryonlybam/ --popmap /cuiwang/nodelete/rawRAD/process_radtags/primaryonlybam/popmap20240123.txt --out-path /cuiwang/nodelete/rawRAD/populations --rm-pcr-duplicates -T 40 -X "populations:--vcf-all" -X "populations:--plink" -X "populations:--treemix" -X "populations:--structure" -X "populations:--radpainter"
#use raxml to build a phylogeny, follow https://rdtarvin.github.io/IBS2019_Genomics-of-Biodiversity/main/2019/08/05/09-raxml-epi.html
cd /cuiwang/nodelete/rawRAD/phylogeny
module load python-data/3.10-23.11
#use the first snp per loci
ln -s /cuiwang/nodelete/rawRAD/populations/firstsnpperloci/populations.snps.vcf .
curl -LO https://raw.githubusercontent.com/edgardomortiz/vcf2phylip/master/vcf2phylip.py
python vcf2phylip.py -i ./populations.snps.vcf
curl -LO https://raw.githubusercontent.com/btmartin721/raxml_ascbias/master/ascbias.py
module load python-data
export PYTHONUSERBASE=/projappl/my-python-env
module load biopythontools

vcftools --vcf populations.snps.vcf --out ref_structure --min-alleles 2 --max-alleles 2 --max-missing 0.5 --maf 0.05 --recode
#remove OG
vcftools --vcf ../ref_structure.recode.vcf --remove-indv Y56 --recode
python vcf2phylip.py -i ./out.recode.vcf
python ascbias.py -p ./populations.min4.phy
raxml-ng --all --msa out.phy --model GTR+ASC_LEWIS --tree pars{10} --bs-trees 100 --threads 8

#chloroplast
cd /cuiwang/nodelete/rawRAD/Chloroplast
module load biokit 
for i in ../process_radtags/fastq/*\.1.fq.gz; do
bwa mem -t 40 Chloroplst_genome.fasta $i ${i%.1.fq.gz}.2.fq.gz | samtools sort -@ 40 -o ${i%.1.fq.gz}.bam
done
#filtering out unmapped reads and supplementary reads
for i in `ls *.bam | sed 's/\.bam//'`; do samtools view -@40 -F 2052 -bo $i".filtered.bam" $i".bam"; done
#call SNPs
cd /cuiwang/nodelete/rawRAD/Chloroplast/filtered
module load stacks
ref_map.pl --samples ./ --popmap /cuiwang/nodelete/rawRAD/process_radtags/primaryonlybam/popmap20240123.txt --out-path ./ --rm-pcr-duplicates -T 40 -X "populations:--vcf" 
vcftools --vcf populations.snps.vcf --out chloroplast --remove-indv Y56 --recode
python vcf2phylip.py -i ./chloroplast.recode.vcf
python ascbias.py -p ./chloroplast.recode.min4.phy
#raxml-ng-mpi --all --msa populations.snps.min4.phy --model GTR+I+G4m --tree pars{10} --bs-trees 100 
raxml-ng --all --msa out.phy --model GTR+ASC_LEWIS --tree pars{10} --bs-trees 100 --threads 8


#cyrus 
cd /data2/Cui/Phragmites/RAD/structure/missing0.5maf05
paste -d'\t' K10.10.meanQ popmap.txt | awk '{ print $11"\t"$12"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10 }' > Popu_K10.txt 
awk '{ print $1"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$11}' Popu_K10.txt >Popu_K10_noempty.txt

module load  r-env
start-r
.libPaths(c("/projappl/project_2001259/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
library(treeio)
library(ggtree)
library(ape)
library(ggrepel)
library(plyr)
library(reshape2)
library(ggstance)
library(ggplotify)
library(wesanderson)
library(viridis)
#library(conStruct)
tree<-read.tree("out.recode.min4.phy.raxml.support")
#tree<-read.tree("/scratch/project_2001259/cuiwang/aligntofinal/RADphylo/RAxML_bipartitionsBranchLabels.out.withinvariant")
to_drop <- "Y56"
OG<-c("Y17_1","E12","Y27")
#OG<-"Y56"
tree_reduced <- drop.tip(tree, to_drop)
options(scipen = 999)#prevent R from converting number to scientific notations
structure<-read.csv("./Popu_K10_noempty.txt", header=TRUE, sep="\t", dec=".")
x<-structure[,c(1:7)]
rooted.tree<-root(tree_reduced,outgroup=OG)
#rooted.tree<-root(tree,outgroup=OG)
admix.props<- as.matrix(structure[c(1:88),c(2:7)], header=T)
xlong <- ddply(melt(x, id.vars = 'Ind'), .(Ind), mutate, prop = value / sum(value))
#xlong <- ddply(melt(x, id.vars = 'V1'), .(V1), mutate, prop = value / sum(value))
calibration = makeChronosCalib(rooted.tree)
#https://github.com/simjoly/CourseComparativeMethods/blob/master/lecture2/PhylogeneticTree.Rmd,https://pedrohbraga.github.io/PhyloCompMethods-in-R-workshop/PhyloCompMethodsMaterial.html,http://phylobotanist.blogspot.com/2018/04/time-calibrated-or-at-least-ultrametric.html
b<-chronos(rooted.tree,model = "relaxed",calibration = calibration,control = chronos.control(dual.iter.max = 100))
#is.ultrametric(b) #[1] TRUE
c<-ladderize(b)
#plot(c, edge.width = 2, cex = 0.5)
#If you want to use your chronogram for other analyses, you might have to convert it into a `phylo` object.https://github.com/simjoly/CourseComparativeMethods/blob/
c <- read.tree(text=write.tree(c))
#bootstrap + geom_nodelab(aes(label=label))
p <- ggtree(c)+ geom_tiplab(size=2, align=TRUE, linesize=.5)
p2<- p + geom_facet(panel = "Admixture", data = xlong, geom = ggstance::geom_barh, 
                aes(x = prop, fill = variable), stat = "identity", width = 1.0 ) + scale_fill_manual(values = c("#7301A8FF","#00336FFF","#CD4071FF","#35B779FF","#FB9A06FF","#F1605DFF")) 
#p2<- p + geom_facet(panel = "Admixture", data = xlong, geom = ggstance::geom_barh, aes(x = value, fill = variable), stat = "identity", width = 1.0 ) + scale_fill_manual(values = c("#CD4071FF", "#7301A8FF","#00336FFF","#F1605DFF","#FB9A06FF","#35B779FF","grey")) 
#https://yulab-smu.top/treedata-book/chapter12.html#facet_widths
facet_widths(p2, widths = c(3, 1))
#Usland#7301A8FF,#AU##F1605DFF,#CN#FB9A06FF,#EU#35B779FF,#MED#00336FFF,#NA#CD4071FF
#adjust the windiows to be bigger to save a nice pic
ggsave("Figure2_BC.pdf", device = pdf())
