# Required Libraries
library(ggplot2)
library(readr)
library(reshape2) # for melt()
library(ggsci)
library(wesanderson)


# Reading the Data
data <- read_csv("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/gradient_analysis/sample_list_iso.csv")
start_col <- which(names(data) == "Adrenal Gland")
end_col <- which(names(data) == "Uterus")

heatmap_data <- as.matrix(data[, start_col:end_col])


# Melting the Data for ggplot
melted_data <- melt(heatmap_data)
melted_data$Var1 <- with(melted_data, factor(Var1, levels = rev(unique(Var1))))


# Calculating the mid-positions for each sample type label
sample_types <- unique(data$sample_type)
num_rows_per_sample <- nrow(data) / length(sample_types)
positions <- (1:length(sample_types)) * num_rows_per_sample - (num_rows_per_sample / 2)

# Plotting the heatmap
pal <- wes_palette("Zissou1", 100, type = "continuous")

p <- ggplot(melted_data, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = pal,limits = c(0, 0.5), name = "Proportion expressed") +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        plot.title = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10))+
  theme(plot.margin = unit(c(1, 1, 1,3), "lines"))

p
ggsave("/Users/zhaosen/Library/CloudStorage/OneDrive-BaylorCollegeofMedicine/Project/Deep-RNAseq/gradient_analysis/iso.pdf",plot = p,
       width = 200,height = 100,units ="mm")


