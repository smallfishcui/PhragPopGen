cd /cuiwang/nodelete/rawRAD/populations/firstsnpperloci
mkdir treemix
#treemix file 
populations -P ../ -O ./treemix -M ./popmap.txt --treemix  --threads 40 --min-maf 0.05 --write-single-snp -R 0.5
awk -F' ' '{ print $1" "$3" "$4" "$5" "$6" "$7 }' populations.treemix > populations_noadmix.treemix
module load treemix
module load phylip/3.697
module load parallel/20210922
#remove the first line, comments
#gzip populations_noadmix.treemix
sh Step1_TreeMix.sh populations_noadmix.treemix.gz 8 100 Usnat 500 consense sixgroup 1 10 10
#disable the evonno method
Rscript Step2and4_TreeMix.R
#optM, choose m=2, two migration events
sh Step3_TreeMix.sh populations_noadmix.treemix.gz 8 100 Usnat 500 2 sixgroup 30 sixgroup_constree.newick consense
Rscript Step2and4_TreeMix.R
