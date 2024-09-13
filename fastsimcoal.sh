#tutorial https://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/fastsimcoal2-activity/
cd /cuiwang/nodelete/rawRAD/populations/firstsnpperloci
grep 'AU' popmap.txt > ./fastsimcoal/AU.txt
grep 'CN' popmap.txt > ./fastsimcoal/CN.txt
grep 'EU' popmap.txt > ./fastsimcoal/EU.txt
grep 'Med' popmap.txt > ./fastsimcoal/Med.txt
grep 'Usland' popmap.txt > ./fastsimcoal/Usland.txt
grep 'Usnat' popmap.txt > ./fastsimcoal/Usnat.txt

./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/AU.txt --preview -a --unfolded #select (20, 34526)
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/AU.txt -a --unfolded --proj 20 -o ./fastsimcoal
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/CN.txt --preview -a --unfolded #choose (16, 33642)
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/CN.txt -a --unfolded --proj 16 -o ./fastsimcoal
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/EU.txt --preview -a --unfolded # choose (66, 55964)
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/EU.txt -a --unfolded --proj 66 -o ./fastsimcoal
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/Med.txt --preview -a --unfolded #choose (22, 33948)
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/Med.txt -a --unfolded --proj 22 -o ./fastsimcoal
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/Usland.txt --preview -a --unfolded #choose (10, 22001)
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/Usland.txt -a --unfolded --proj 10 -o ./fastsimcoal
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/Usnat.txt --preview -a --unfolded #choose ((6, 8403))
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./fastsimcoal/Usnat.txt -a --unfolded --proj 6 -o ./fastsimcoal

#stairway plot
java -cp stairway_plot_es Stairbuilder AU.blueprint
java -cp stairway_plot_es Stairbuilder CN.blueprint
java -cp stairway_plot_es Stairbuilder EU.blueprint
java -cp stairway_plot_es Stairbuilder Med.blueprint
java -cp stairway_plot_es Stairbuilder Usland.blueprint
java -cp stairway_plot_es Stairbuilder Usnat.blueprint

#plot
awk '{ print "AU",$6/1000,$7/1000,$8/1000,$9/1000 }' ./AU/AU.final.summary > AU.final.txt
awk '{ print "CN",$6/1000,$7/1000,$8/1000,$9/1000 }' ./CN/CN.final.summary > CN.final.txt
awk '{ print "EU",$6/1000,$7/1000,$8/1000,$9/1000 }' ./EU/EU.final.summary > EU.final.txt
awk '{ print "Med",$6/1000,$7/1000,$8/1000,$9/1000 }' ./Med/Med.final.summary > Med.final.txt
awk '{ print "Usland",$6/1000,$7/1000,$8/1000,$9/1000 }' ./Usland/Usland.final.summary > Usland.final.txt
awk '{ print "Usnat",$6/1000,$7/1000,$8/1000,$9/1000 }' ./Usnat/Usnat.final.summary > Usnat.final.txt
#remove the first lines in all the files
cat AU.final.txt CN.final.txt EU.final.txt Med.final.txt Usland.final.txt Usnat.final.txt > stairway.plot.txt
#add Group year(1k) Ne_median(1k) Ne_2.5%(1k) Ne_97.5%(1k) to the header

module load  r-env
start-r
.libPaths(c("/projappl/project_2001259/project_rpackages", .libPaths()))
libpath <- .libPaths()[1]
library(survival)
library(broom)
library(ggplot2)
library(utile.visuals)
library(brms)
a<-read.table("stairway.plot.txt", sep=" ", header=T)
#awk '{print $1,2.5*log($2)/log(10),2.5*log($3)/log(10) }' stairway.plot.txt > stairway.plot.log10.txt
#nano the header
#ggplot+log10 scale
a<-read.table("stairway.plot.txt", sep=" ", header=T)
pdf("stairway.plot.log10.pdf")
ggplot(
  data = a,
  mapping = aes(x = year.1k., y = Ne_median.1k.)
) +
  geom_step(aes(color = Group), size=3) +
  labs(x = 'log10year(1k)', y = 'log10Ne(1k)') + scale_x_log10(n.breaks =10) +
    scale_y_log10()+ 
    scale_color_manual(
    values = c("AU"="#F1605DFF","Usnat"="#CD4071FF","Med"="#00336FFF","EU"="#35B779FF","CN"="#FB9A06FF","Usland"="#7301A8FF"),
    aesthetics = c('colour', 'fill')) +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.text=element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold")) + theme(legend.title = element_text(size=14,face="bold"), legend.text = element_text(size=14,face="bold"))

dev.off()
ggsave("stairway.plot.log10.pdf")

#fastsimcoal
grep -v 'admix' popmap.txt > allpurepop.txt
./easySFS/easySFS.py -i ./Ancestral/output.AA.swap.vcf.gz -p ./allpurepop.txt -a --unfolded --proj 20,16,66,22,10,6 -o ./fastsimcoal
cd 
mkdir AU_CN
PREFIX=USland_USnat
mv ${PREFIX}.* AU_CN
cd AU_CN
#PREFIX=AU_CN
for i in {1..100}
 do
   mkdir run$i
   cp ${PREFIX}.tpl ${PREFIX}.est ../fsc28_linux64/fsc28 ${PREFIX}_jointDAFpop1_0.obs run$i"/"
   cd run$i
   ./fsc28 -t ${PREFIX}.tpl -e ${PREFIX}.est -n 200000 -M -d -L 40 -c40
   cd ..
 done
#cat run{1..100}/${PREFIX}/${PREFIX}.bestlhoods | grep -v MaxObsLhood | awk '{print NR,$6}' | sort -k 2 | head, not sure if this is right
bash fsc-selectbestrun.sh

#use all the SNP sites to include invariant sites
#cd /scratch/project_2001259/cuiwang/nodelete/rawRAD/populations/all
#./easySFS/easySFS.py -i ./output.AA.swap.vcf.gz -p ./allpurepop.txt -a --unfolded --proj 20,16,66,22,10,6 -o ./fastsimcoal -f
#on mahti

##################final#########################
#In total, there are 34997178 variant sites, select first snp per loci, we get 643615 sites. (54.37x)
#these represent a length of 911125290 sites, so the first snp per loci represent 16757868 sites.
ln -s /cuiwang/nodelete/rawRAD/populations/firstsnpperloci/Ancestral/output.AA.swap.vcf.gz .
./easySFS/easySFS.py -i ./output.AA.swap.vcf.gz -p ./allpurepop.txt -a --unfolded --proj 20,16,66,22,10,6 -o ./fastsimcoal_allpop --total-length 16757868
./easySFS/easySFS.py -i ./output.AA.swap.vcf.gz -p ./Med_EU_CN.txt -a --unfolded  --preview -o ./fastsimcoal --total-length 16757868
./easySFS/easySFS.py -i ./output.AA.swap.vcf.gz -p ./Med_EU_CN.txt -a --unfolded --proj 12,46,14 -o ./fastsimcoal --total-length 16757868 # sites 53952,101864,52665
mv fastsimcoal fastsimcoalMed_EU_CN
./easySFS/easySFS.py -i ./output.AA.swap.vcf.gz -p ./CN_EU_Med_USnat.txt -a --unfolded --preview -o ./fastsimcoal --total-length 16757868 # sites 53952,101864,52665
./easySFS/easySFS.py -i ./output.AA.swap.vcf.gz -p ./CN_EU_Med_USnat.txt -a --unfolded --proj 12,46,14,6 -o ./fastsimcoal --total-length 16757868 # sites 53952,101864,52665,8403


#calculate AIC
AIC=2*k-2*(bestlhoods$MaxEstLhood/log10(exp(1)))
#found the best model is current gene flow between both eu and med ,and eu and cn

#after the initial runs, run the best models for 100 times
for i in {1..100}
do
 ./fsc28 -i ${PREFIX}_maxL.par -n1000000 -d 
 # Fastsimcoal will generate a new folder called ${model}_maxL and write files in there
 sed -n '2,3p' ${PREFIX}_maxL/${PREFIX}_maxL.lhoods  >> ${PREFIX}.lhoods
 rm -r ${PREFIX}_maxL/
done


