setwd()
kegg_up <- read.csv("KEGG_up.csv")
library(ggplot2)
plot_up <- ggplot(kegg_up,aes(y=Description,x=Count))+
  geom_point(aes(size=Count,color=-1*PValue))+#      
  scale_color_gradient(low="green",high = "red")+
  labs(color=expression(PValue,size="Count"), 
       x="Gene Number",y="Pathways",title="KEGG Pathway Enrichment for Upregulated Genes")+
  theme_bw()
ggsave("KEGG_up.png", plot_up, height = 5, width = 4 * 2)
library(gridExtra)
gridExtra::grid.arrange(plot1, plot2, ncol=2)
