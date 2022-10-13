
# Script 6 ----------------------------------------------------------------
# Figure 1A and 5A

# Skin Media Manuscript
# Phylogenetic tree and tanglegram
# Written by Mary Hannah Swaney
# University of Wisconsin - Madison

# Load required packages and data -----------------------------------------

library(ggtree)
library(ape)
library(dendextend)

tree <- read.tree("/Users/mhswaney/Lab_data_permanent/Skin_media_project/SkinStrains_PhylogeneticTree/SkinMedia16S.tree")

# 16S phylogeny -----------------------------------------------------------

# root the tree
rooted <- root(tree, which(tree$tip.label == "Halobacterium_salinarum")) %>% drop.tip(c("Halobacterium_salinarum"))

# re-label tips
rooted$tip.label[rooted$tip.label=="Kocuria_rhizophila_LK221"] <- "Kocuria rhizophila LK221"
rooted$tip.label[rooted$tip.label=="Kocuria_marina_A_LK478"] <- "Kocuria marina_A LK478"
rooted$tip.label[rooted$tip.label=="Dermacoccus_nishinomiyaensis_LK1128"] <- "Dermacoccus nishinomiyaensis LK1128"
rooted$tip.label[rooted$tip.label=="Dermabacter_hominis_LK522"] <- "Dermabacter hominis LK522"
rooted$tip.label[rooted$tip.label=="Micrococcus_luteus_LK410"] <- "Micrococcus luteus LK410"
rooted$tip.label[rooted$tip.label=="Corynebacterium_kefirresidentii_LK1134"] <- "Corynebacterium kefirresidentii LK1134"
rooted$tip.label[rooted$tip.label=="Corynebacterium_glucuronolyticum_LK488"] <- "Corynebacterium glucuronolyticum LK488"
rooted$tip.label[rooted$tip.label=="Corynebacterium_amycolatum_LK19"] <- "Corynebacterium amycolatum LK19"
rooted$tip.label[rooted$tip.label=="Dietzia_cinnamea_LK439"] <- "Dietzia cinnamea LK439"
rooted$tip.label[rooted$tip.label=="Microbacterium_sp._LK369"] <- "Microbacterium sp. LK369"
rooted$tip.label[rooted$tip.label=="Staphylococcus_epidermidis_LK717"] <- "Staphylococcus epidermidis LK717"
rooted$tip.label[rooted$tip.label=="Staphylococcus_epidermidis_LK593"] <- "Staphylococcus epidermidis LK593"
rooted$tip.label[rooted$tip.label=="Staphylococcus_aureus_USA300"] <- "Staphylococcus aureus USA300"
rooted$tip.label[rooted$tip.label=="Citrobacteri_freundii_LK704"] <- "Citrobacter freundii LK704"
rooted$tip.label[rooted$tip.label=="Klebsiella_pneumoniae_LK469"] <- "Klebsiella pneumoniae LK469"
rooted$tip.label[rooted$tip.label=="Sphingobacterium_hotanense_LK485"] <- "Sphingobacterium hotanense LK485"

final <- ggtree(rooted) + geom_tiplab(size=3, align=TRUE) +
  xlim(0,1.5)

# Figure 1A
final

#ggsave(plot = final, filename =  "/Users/mhswaney/Kalan_Lab/Manuscripts/SkinMedia/Figures/Figure1/16STree_SkinMediaStrains.pdf",
      # device = "pdf",width = 4, height=6,units = "in")

# Tanglegram --------------------------------------------------------------

# Create tanglegram
dendrogram <- chronos(rooted)
tanglegram(hclust_avg, dendrogram, sort=TRUE, # hclust_avg is from script 5
           common_subtrees_color_branches = TRUE,
           highlight_branches_col = TRUE)

# create list of the two dendrograms
both <- dendlist(as.dendrogram(dendrogram), as.dendrogram(hclust_avg))

both %>% entanglement # calculate entanglement value

# plot tanglegram - Figure 5A
both %>% dendextend::untangle(method = "step2side") %>%
  tanglegram( common_subtrees_color_lines = TRUE,
              common_subtrees_color_branches = TRUE,
              highlight_branches_col = FALSE,
              highlight_branches_lwd = FALSE,
              highlight_distinct_edges = FALSE,
              lab.cex = 1.2, margin_inner = 22)

