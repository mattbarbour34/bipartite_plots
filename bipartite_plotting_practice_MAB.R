library(bipartite)
library(GGally)
source("~/Documents/bipartite_plots/functions/bip_briatte.R")
library(GGally)
library(dplyr)
library(gridExtra)

genotype_gall_parasitoid_network <- read.csv('~/Documents/Genotype_Networks/data/genotype_gall_parasitoid_network_data.csv')
rownames(genotype_gall_parasitoid_network) <- genotype_gall_parasitoid_network$X 
genotype_gall_parasitoid_network <- select(genotype_gall_parasitoid_network, -X)

plotweb(genotype_gall_parasitoid_network)

web <- genotype_gall_parasitoid_network
ca <- cca(web)
web <- web[order(summary(ca)$sites[, 1], decreasing = TRUE), 
           order(summary(ca)$species[, 1], decreasing = TRUE)]
low.order <- order(summary(ca)$sites[, 1], decreasing = TRUE)
high.order <- order(summary(ca)$species[, 1], 
                    decreasing = TRUE)

module.info <- read.csv('~/Documents/Genotype_Networks/data/module_info_genotype_gall_parasitoid_network.csv')

net = bipartite.network(genotype_gall_parasitoid_network, modes = c("Genotypes", "Gall-Parasitoid Interactions"), directed = FALSE)
set.vertex.attribute(net, "modules", c(module.info$Module.ID[1:25], module.info$Module.ID[26:39]))

layout.bipartite <- matrix(nrow = 39, ncol = 2)
layout.bipartite[ ,1] <- c(low.order, high.order)

x = 5
layout.bipartite[ ,2] <- c(rep(1, 25), rep(1*x, 14))

layouts <- plot(net, coord = layout.bipartite)


ggnet(net,
      coord = layout.bipartite,
      segment.size = edge.weights(genotype_gall_parasitoid_network, 10), 
      segment.alpha = .5,
      label = TRUE, color = "black",
      node.group = get.vertex.attribute(net, "mode"))
# Don't think this is a promising visualization....

net_melt <- mutate(genotype_gall_parasitoid_network, Genotype = factor(rownames(genotype_gall_parasitoid_network), levels = rev(c("O","A","Z","H","J","S","K","N","Y","G","D","V","X","P","B","W","R","Q","E","T","F","M","L","I","*")), ordered = TRUE)) 

net_melt <- melt(net_melt) 
  
net_melt <- mutate(net_melt, gall.ptoid.interaction = factor(variable, levels = c("vLG_Tory","rG_Tory","vLG_Mymarid","vLG_Eulo","SG_Platy","vLG_Platy","rG_Mesopol","rG_Platy","rG_Lestodip","rG_Eulo","rsLG_Eury","aSG_Tory","vLG_Mesopol","rsLG_Lathro"), ordered = TRUE))

theme_bipartite <- theme(axis.text = element_text(size=12),
                         axis.text.x = element_text(angle = 45,
                                                    hjust = 1),
                         axis.ticks = element_blank(),
                         panel.background = element_blank())

module.color = "steelblue"
ggplot(data = net_melt,
       aes(x = gall.ptoid.interaction, y = Genotype, fill = value)) + 
  geom_tile(color = "black") + 
  #geom_text(label = net_melt$value) +
  scale_fill_gradient(low = "white",
                      high = "red") +
  xlab("Gall-Parasitoid Interaction") +
  theme_bipartite + 
  annotate("rect", xmin=0.5, xmax=4.5, ymin=17.5, ymax=Inf, 
           color=module.color, alpha = 0, size = 2) +
  annotate("rect", xmin=4.5, xmax=5.5, ymin=16.5, ymax=17.5, 
           color=module.color, alpha = 0, size = 2) + 
  annotate("rect", xmin=5.5, xmax=7.5, ymin=10.5, ymax=16.5, 
           color=module.color, alpha = 0, size = 2) + 
  annotate("rect", xmin=7.5, xmax=10.5, ymin=7.5, ymax=10.5, 
           color=module.color, alpha = 0, size = 2) + 
  annotate("rect", xmin=10.5, xmax=14.5, ymin=0.5, ymax=7.5, 
           color=module.color, alpha = 0, size = 2)