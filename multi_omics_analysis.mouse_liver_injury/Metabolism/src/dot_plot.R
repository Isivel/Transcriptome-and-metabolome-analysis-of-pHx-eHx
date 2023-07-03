
library(ggplot2)
library(enrichR)
library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)
library(stringr)
library(reshape2)

##### GSEA
dotplot_ylab <- function(x,showCategory = 20, wi = 50, term_cex = 8, term = NA){
  x.df <- as.data.frame(x)
  
  if(nrow(x.df)>showCategory){
    x.df <- x.df[1:showCategory,]
  }
  x.df$GR <- as.numeric(gsub("^(\\d+)\\/\\d+$", "\\1",x.df$GeneRatio)) / as.numeric(gsub("^\\d+\\/(\\d+)$", "\\1",x.df$GeneRatio))
  #x.df$oddRatios <-  unlist(lapply(as.list(x.df$oddRatios), function(x) eval(parse(text=x))))
  x.df$Description <- factor(x.df$Description,levels = x.df$Description[order(x.df$GR,decreasing = F)])
  
  #x.df <- x.df[order(x.df$Count,decreasing = F),]
  p<-ggplot(x.df, aes(x=GR, y=Description,color=pvalue)) +
    geom_point(aes(size = Count))+scale_color_gradient(low="red", high="blue") +
    labs(x = "GeneRatio", title = term) +
    scale_y_discrete(labels=function(x) str_wrap(x,width=wi))+
    theme(axis.text.y = element_text( size = term_cex),
          plot.title = element_text(hjust = .5)) +
    guides(color = guide_legend(order = 0),
           size = guide_legend(order = 1))
}


