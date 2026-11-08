library(readxl)
library(tidyverse)

data=read_xlsx('Suppl table.xlsx')%>%pivot_longer(col=-1,names_to = 'Predict',values_to = 'Value')
colnames(data)[1] = 'Truth'
head(data)


data$Truth=factor(data$Truth,levels = c('B cells','CD14 Mono','CD16 Mono','CD4 Tcells','CD8 Tcells','NK cells','pDC','DC'))
data$Predict=factor(data$Predict,levels = rev(c('B cells','CD14 Mono','CD16 Mono','CD4 Tcells','CD8 Tcells','NK cells','pDC','DC','cDC','CMP','Basophils','Eosinophils','Erythrocytes','GMP','HSC','Macrophages','Megakaryocytes','MEP','MPP','Neutrophils','Platelets')))



SupFig2A=ggplot(data, aes(y = Predict, x = Truth, fill = Value)) +
  geom_tile(color = "black") +  
  geom_text(aes(label = round(Value, 2)),size=2) +
 scale_fill_gradient(high = "#d73027", low = "#FEFAFA")+
  theme_minimal() +
  labs(title = NULL, x = NULL, y = NULL) +
  theme(
      axis.title.x = element_blank(),  
      axis.text.x = element_blank(),   
      legend.title = element_blank()  
  )
SupFig2A
ggsave(file='SupFig2A.pdf',P,width = 4,height = 5.5)