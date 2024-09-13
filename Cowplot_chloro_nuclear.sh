cd /cuiwang/nodelete3/rawRAD/figures
ln -s /cuiwang/nodelete3/rawRAD/phylogeny/out.recode.min4.phy.raxml.support .
ln -s /cuiwang/nodelete3/rawRAD/Chloroplast/filtered/out.phy.raxml.support .
mv out.recode.min4.phy.raxml.support nuclear.tre
mv out.phy.raxml.support chloroplast.tre

module load r-env
start-r
library(dplyr)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(ape)
library(phytools)
library("ape")

# Read the trees from files
tree1 <- read.tree("nuclear.tre", tree.names = NULL)
tree2 <- read.tree("chloroplast.tre", tree.names = NULL)

# Define the outgroup for rooting both trees
outgroup_cluster <- c("Y27", "E12", "Y17_1")

# Root both trees using the defined outgroup
rooted_tree1 <- root(tree1, outgroup = outgroup_cluster)
rooted_tree2 <- root(tree2, outgroup = outgroup_cluster)
assoc=cbind(rooted_tree1$tip.label, rooted_tree1$tip.label)

pdf("face_to_face_trees_with_links.pdf", width = 10, height = 15)  # Output to PDF, adjust width and height as needed
cophyloplot(rooted_tree1, rooted_tree2,assoc =assoc, space = 50, gap = 10,length.line = 0.5, lwd = 2, col = "blue")  # Adjust space as needed
dev.off()  # Close the PDF device

