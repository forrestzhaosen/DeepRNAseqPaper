library(readr)
library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(ggsci)
library(gridExtra)

setwd("/Users/zhaosen/PUMCH Dropbox/赵森/Baylor/Project/Deep-RNAseq/spliceVault_construction")
data <- read_tsv("dot_data.txt")

top1 <- data %>%
  filter(gold_rank==1)%>%
  mutate(rank_class = case_when(
    rank == 1 ~ "top1",
    rank %in% 2:3 ~ "top3",
    is.na(rank) ~ "missing",
    TRUE ~ "other"  # Handles any other cases if needed
  )) %>%
  mutate(rank_class = factor(rank_class, levels = c("missing","top3","top1"))) %>%
  group_by(class,rank_class)%>%
  summarise(rank_proportion = n())

p1<-ggplot(top1, aes(fill=rank_class, y=rank_proportion, x=class)) + 
  geom_bar(position="fill", stat="identity") +
  labs(x=NULL,y="Proportion")+
  scale_fill_nejm(name = NULL)+
  theme_classic() +
  theme(text = element_text(size = 20))
p1

multi <- data %>%
  filter(multi==1)


# Creating the dot plot with primary grouping by Gene and secondary grouping by Class
p2 <- ggplot(multi, aes(x = class, y = rank, color = class)) +
  geom_point(position = position_dodge(width = 0.1)) +  # Dodge position to handle overlapping points
  facet_wrap(~ Gene, scales = "free_x", ncol = 3, labeller = label_both) +  # Two columns, label both gene and class
  scale_color_nejm(name = NULL)+
  labs(title = NULL,
       x = NULL,
       y = "Rank") +
  theme_classic() +
  theme(text = element_text(size = 20),axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x labels for readability
p2

p=do.call(grid.arrange, c(list(p1,p2), ncol=2))
ggsave("splicevault.pdf",plot =p,
       width = 360,height = 130,units ="mm")
