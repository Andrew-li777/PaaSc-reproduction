library(tidyverse)
library(ggpubr)
data=read.table('/sibcb1/bioinformatics/fengyan/Project/fig3/Bcell/Marker/Bcell_marker_auc.txt')%>%pivot_longer(.,cols=1:8,names_to ='Method',values_to='AUC')%>%
  group_by(Method) %>%
  mutate(median_auc = median(AUC, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Method = factor(Method, levels = unique(Method[order(-median_auc)])))
ggplot(data, aes(x = Method, y = AUC, color = Method)) +
          geom_boxplot(outlier.shape = NA) +
          geom_jitter(size=5)+
          guides(fill = "none", color = "none")+
         scale_color_manual(values = c('steelblue','black','#56A902','#2e409a','#fcbe32','#e97f02','#f15c5c','#fee3ce'))+
          labs( x = NULL, y = 'Score',title = NULL) +
          theme_bw() +
          theme(
                axis.ticks.y = element_line(color = "black",size=15), 
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,, size = 15)#,
             
          )